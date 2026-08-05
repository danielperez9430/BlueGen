#!/usr/bin/env python3
"""
PGS Population Calibration — scores all 1000G samples with each PGS score,
builds per-population distributions, and calibrates target sample.

Usage:
  python3 pgs_population_calibrate.py \
    --bfile reference/1000G_full/1000G_full \
    --pop-panel reference/1000G_full/population_panel.txt \
    --pgs-dir prs/pgs_scores \
    --sample-prs prs/pgs_scores/pgs_results.csv \
    --output-dir prs/pgs_scores/
"""

import sys
import os
import subprocess
import json
from datetime import datetime
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import scipy.stats

MAX_RELIABLE_SNPS = 500_000  # matches the platform's documented reliability cutoff


def read_score_metadata(pgs_dir: Path, pgs_id: str) -> dict:
    """Read trait name + variant count from a PGS score file's own header
    comments (#trait_reported=...) and body, so the calibration report can
    show a human-readable trait name instead of just the PGS ID."""
    score_dir = Path(pgs_dir) / pgs_id
    candidates = sorted(score_dir.glob(f"{pgs_id}*.txt")) if score_dir.exists() else []
    if not candidates:
        return {"trait": pgs_id, "n_snps": 0}
    raw_file = candidates[0]
    trait = pgs_id
    n_snps = 0
    header_done = False
    with open(raw_file, errors="replace") as fh:
        for line in fh:
            if line.startswith("#trait_reported="):
                trait = line.split("=", 1)[1].strip() or pgs_id
            elif line.startswith("#"):
                continue
            elif not header_done:
                header_done = True  # column header row
            else:
                n_snps += 1
    return {"trait": trait, "n_snps": n_snps}

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bfile", required=True)
    p.add_argument("--pop-panel", required=True)
    p.add_argument("--pgs-dir", required=True)
    p.add_argument("--sample-prs", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plink", default="plink")
    p.add_argument("--threads", default="8")
    args = p.parse_args()

    bfile = Path(args.bfile)
    pgs_dir = Path(args.pgs_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plink = Path(args.plink)
    ref_dir = out_dir / "ref_distributions"
    ref_dir.mkdir(exist_ok=True)

    # Load population panel
    pop = pd.read_csv(args.pop_panel, sep="\t", dtype=str)
    pop_map = dict(zip(pop["sample"], pop["super_pop"]))
    print(f"Population panel: {len(pop)} samples across {len(pop['super_pop'].unique())} super-pops")

    # Load sample PRS
    sample_df = pd.read_csv(args.sample_prs, dtype=str)
    sample_df["prs_raw"] = sample_df["prs_raw"].astype(float)
    print(f"Sample PRS: {len(sample_df)} entries ({sample_df['pgs_id'].nunique()} scores)")

    # Per-PGS calibration
    results = []
    pgs_files = sorted(Path(pgs_dir).glob("PGS*/PGS*_clean.score"))

    # PLINK reloads and re-indexes the full genome-wide reference (84M
    # variants) on every invocation - doing that once per PGS score (dozens
    # of them) made a full run take 1.5-2+ hours. Extract the union of
    # variant IDs actually needed across every score ONCE into a much
    # smaller working bfile, then --score against that instead of the raw
    # reference on each iteration.
    needed_ids = set()
    for score_file in pgs_files:
        with open(score_file) as fh:
            for line in fh:
                vid = line.split("\t", 1)[0].strip()
                if vid:
                    needed_ids.add(vid)
    print(f"\n{len(needed_ids)} unique variant IDs needed across {len(pgs_files)} scores")

    score_bfile = bfile
    if needed_ids:
        extract_list = out_dir / "_needed_variants.txt"
        extract_list.write_text("\n".join(sorted(needed_ids)))
        subset_prefix = out_dir / "_reference_subset"
        r = subprocess.run([
            str(plink), "--bfile", str(bfile), "--extract", str(extract_list),
            "--make-bed", "--out", str(subset_prefix),
            "--allow-extra-chr", "--threads", args.threads, "--memory", "16000"
        ], capture_output=True, text=True, timeout=1800)
        if Path(str(subset_prefix) + ".bed").exists():
            score_bfile = subset_prefix
            print(f"  Reference subset built: {subset_prefix}")
        else:
            print(f"  Reference subset extraction failed ({r.stderr[-200:]!r}), "
                  "falling back to the full reference per-score (slow)")

    for i, score_file in enumerate(pgs_files):
        pgs_id = score_file.parent.name
        print(f"\n[{i+1}/{len(pgs_files)}] {pgs_id} ({pgs_id})")

        out_prefix = ref_dir / pgs_id

        # Run PLINK --score against 1000G (or the pre-extracted subset above)
        r = subprocess.run([
            str(plink), "--bfile", str(score_bfile),
            "--score", str(score_file), "1", "2", "3",
            "--out", str(out_prefix),
            "--allow-extra-chr", "--threads", args.threads, "--memory", "16000"
        ], capture_output=True, text=True, timeout=600)

        prof_path = Path(str(out_prefix) + ".profile")
        if not prof_path.exists():
            print(f"  ❌ No profile — skipping")
            continue

        # Parse 1000G profile
        ref_scores = pd.read_csv(prof_path, sep=r"\s+", dtype={"IID": str})
        ref_scores["super_pop"] = ref_scores["IID"].map(pop_map)
        ref_scores["SCORE"] = ref_scores.get("SCORE", ref_scores.get("SCORESUM", 0)).astype(float)

        # Build per-population distributions
        dist = {}
        for sp in ["EUR", "AFR", "EAS", "SAS", "AMR"]:
            pop_scores = ref_scores[ref_scores["super_pop"] == sp]["SCORE"]
            if len(pop_scores) < 10:
                continue
            dist[sp] = {
                "n": int(len(pop_scores)),
                "mean": float(pop_scores.mean()),
                "std": float(pop_scores.std()),
                "median": float(pop_scores.median()),
                "p5": float(np.percentile(pop_scores, 5)),
                "p10": float(np.percentile(pop_scores, 10)),
                "p25": float(np.percentile(pop_scores, 25)),
                "p75": float(np.percentile(pop_scores, 75)),
                "p90": float(np.percentile(pop_scores, 90)),
                "p95": float(np.percentile(pop_scores, 95)),
                "min": float(pop_scores.min()),
                "max": float(pop_scores.max()),
            }

        # Save distribution
        with open(ref_dir / f"{pgs_id}_dist.json", "w") as f:
            json.dump(dist, f, indent=2)

        # Calibrate sample
        sample_row = sample_df[sample_df["pgs_id"] == pgs_id]
        if len(sample_row) > 0:
            sample_score = float(sample_row.iloc[0]["prs_raw"])
            # Use EUR as default (most common for our sample)
            eur = dist.get("EUR", {})
            if eur:
                z = (sample_score - eur["mean"]) / eur["std"] if eur["std"] > 0 else 0
                # Previously gated on `'scipy' in dir()`, which only checks
                # local names inside this function - scipy.stats was only
                # ever imported under `if __name__ == "__main__":` at module
                # scope, so that check was always False and percentile was
                # always the hardcoded fallback of 50, regardless of z.
                pctl = scipy.stats.norm.cdf(z) * 100
                risk = "HIGH" if z > 2 else ("ELEVATED" if z > 1 else ("AVERAGE" if abs(z) <= 1 else ("LOW" if z < -1 else "PROTECTIVE")))
            else:
                z, pctl, risk = 0, 50, "UNKNOWN"

            meta = read_score_metadata(pgs_dir, pgs_id)
            results.append({
                "pgs_id": pgs_id,
                "trait": meta["trait"],
                "n_snps": meta["n_snps"],
                "sample_score": sample_score,
                "eur_mean": eur.get("mean", np.nan),
                "eur_std": eur.get("std", np.nan),
                "z_score": round(z, 3),
                "percentile": round(pctl, 1),
                "risk_category": risk,
                "reliable": meta["n_snps"] <= MAX_RELIABLE_SNPS,
                "n_populations": len(dist),
            })

        # Cleanup
        prof_path.unlink(missing_ok=True)
        for ext in [".log", ".nosex"]:
            Path(str(out_prefix) + ext).unlink(missing_ok=True)

    # Save calibrated results
    if results:
        cal = pd.DataFrame(results).sort_values("z_score", ascending=False)
        cal.to_csv(out_dir / "pgs_calibrated.csv", index=False)

        print(f"\n{'='*80}")
        print(f"{'PGS ID':<14} {'Score':>12} {'Z-Score':>10} {'Pctl':>8} {'Risk':>14}  {'EUR mean ± std':>25}")
        print("-" * 80)
        for _, r in cal.iterrows():
            label = "🔴" if r["risk_category"] == "HIGH" else ("🟠" if r["risk_category"] == "ELEVATED" else ("🟡" if r["risk_category"] == "AVERAGE" else "🟢"))
            print(f"{r['pgs_id']:<14} {r['sample_score']:>12.6f} {r['z_score']:>10.2f} {r['percentile']:>7.1f}% {label} {r['risk_category']:<10}  μ={r['eur_mean']:.4f} σ={r['eur_std']:.4f}")
        print(f"\n✅ {len(results)} PGS scores calibrated against 1000G EUR population")
        print(f"📁 Results: {out_dir / 'pgs_calibrated.csv'}")

        # comprehensive_report.py's PGS Catalog section reads this structured
        # report, not the flat CSV above - previously nothing wrote it at all
        # (prs/pgs_scores/pgs_calibration_report.json sat frozen since
        # 2026-07-16), so the report's PGS section silently never reflected a
        # live run either, independent of the CSV going stale.
        report = {
            "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
            "methodology": {
                "reference_panel": f"1000 Genomes Phase 3 ({len(pop)} samples)",
                "populations": ["EUR", "AFR", "EAS", "SAS", "AMR"],
                "calibration_method": "Population-stratified z-score normalization",
                "note": f"Scores with >{MAX_RELIABLE_SNPS:,} SNPs flagged as potentially "
                        "unreliable due to distribution compression",
            },
            "summary": {
                "total_scores": len(results),
                "reliable_scores": int(cal["reliable"].sum()),
                "high_risk": int((cal["risk_category"] == "HIGH").sum()),
                "elevated_risk": int((cal["risk_category"] == "ELEVATED").sum()),
                "low_risk": int((cal["risk_category"] == "LOW").sum()),
            },
            "high_risk_traits": cal[cal["risk_category"] == "HIGH"].to_dict("records"),
            "elevated_risk_traits": cal[cal["risk_category"] == "ELEVATED"].to_dict("records"),
            "low_risk_traits": cal[cal["risk_category"] == "LOW"].to_dict("records"),
            "all_entries": cal.to_dict("records"),
            "reference_distributions_path": str(ref_dir) + "/",
        }
        with open(out_dir / "pgs_calibration_report.json", "w") as f:
            json.dump(report, f, indent=2, default=float)
        print(f"📁 Report: {out_dir / 'pgs_calibration_report.json'}")

if __name__ == "__main__":
    main()
