#!/usr/bin/env python3
"""
PGS population calibration — JOINT scoring of the user and the 1000 Genomes
reference on one identical variant set, then z-score / percentile of the
user against the super-population inferred by the ancestry stage.

Why joint scoring (RELEASE_PLAN 3.0.1, 2026-09-14)
--------------------------------------------------
The previous version scored the user alone on qc/qc_filtered and the
reference on every variant of each score. A WGS-derived PLINK dataset only
contains sites where the sample carries at least one ALT allele (hom-ref
"RefCall" lines are dropped in Stage A), so the user was summed over a
biased subset while 2,504 reference samples were summed over the full set.
The two numbers were not comparable and produced z-scores such as -131
(T2D) or -38 (BMI), with 12 of 52 scores at percentile exactly 0 or 100.

This version:
  1. extracts the score variants from the reference AND from the user,
     restricted to chromosomes the user's data actually covers,
  2. merges both into one dataset (reference first, --keep-allele-order, so
     A2 stays the hg19 REF allele of the 1000G conversion),
  3. --fill-missing-a2: every variant the user lacks becomes homozygous REF
     — the same assumption the curated PRS path makes in 06_prs_compute.py
     (absent from a WGS VCF == homozygous reference),
  4. runs a single `plink --score ... sum` per PGS on the joint dataset, so
     the user and every reference sample are summed over exactly the same
     variants with the same allele matching,
  5. builds per-super-population distributions and calibrates the user
     against `assigned_population` from science/ANCESTRY_MODEL.json
     (EUR only as an explicit, flagged fallback).

Coverage is reported honestly: n_snps_matched / n_snps of the score; a score
is `reliable` only when it is ≤ MAX_RELIABLE_SNPS AND ≥ MIN_COVERAGE of its
variants were in the joint set.

Usage:
  python3 pgs_population_calibrate.py \
    --bfile reference/1000G_full/1000G_full \
    --pop-panel reference/1000G_full/population_panel.txt \
    --pgs-dir pgs \
    --user-bfile qc/qc_filtered \
    --ancestry-json science/ANCESTRY_MODEL.json \
    --output-dir prs/pgs_scores
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

SUPER_POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
MAX_RELIABLE_SNPS = 500_000   # platform's documented reliability cutoff
MIN_COVERAGE = 0.80           # fraction of a score's variants that must be in the joint set
MIN_REF_PER_POP = 10          # minimum reference samples to build a distribution
AUTOSOMES = [str(c) for c in range(1, 23)]


# ── helpers ──────────────────────────────────────────────────────────────────

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


from bluegen.joint import (JointDatasetError, build_joint_dataset, read_fam_ids,  # noqa: E402
                           remove_joint_dataset, run_plink, user_covered_chromosomes)


def load_assigned_population(ancestry_json) -> tuple:
    """Returns (super_pop, source). source is 'inferred' or 'fallback'."""
    if ancestry_json:
        p = Path(ancestry_json)
        if p.exists():
            try:
                with open(p) as fh:
                    pop = str(json.load(fh).get("assigned_population", "")).upper().strip()
                if pop in SUPER_POPS:
                    return pop, "inferred"
                print(f"⚠️  {p}: assigned_population={pop!r} is not a 1000G super-population; falling back to EUR")
            except (OSError, ValueError) as e:
                print(f"⚠️  Could not read {p}: {e}; falling back to EUR")
        else:
            print(f"⚠️  Ancestry file not found: {p}; falling back to EUR")
    else:
        print("⚠️  No --ancestry-json given; falling back to EUR")
    return "EUR", "fallback"


def dedup_score_file(score_file: Path, out_path: Path) -> tuple:
    """PLINK aborts on a duplicated variant ID in a --score file. Clean score
    files use chr:pos as the ID, so multi-allelic sites (two rows, same
    position, different alleles) collide. Drop every row of a duplicated ID
    (PLINK could not disambiguate them by position anyway). Returns
    (path_to_use, n_rows_dropped)."""
    counts = {}
    rows = []
    with open(score_file) as fh:
        for line in fh:
            vid = line.split("\t", 1)[0].strip()
            if not vid:
                continue
            counts[vid] = counts.get(vid, 0) + 1
            rows.append((vid, line))
    dups = {v for v, c in counts.items() if c > 1}
    if not dups:
        return score_file, 0
    kept = [line for vid, line in rows if vid not in dups]
    out_path.write_text("".join(kept))
    return out_path, len(rows) - len(kept)


def score_joint(plink, joint_prefix, score_file, out_prefix, threads, memory) -> pd.DataFrame:
    """One `--score ... sum` over the joint dataset. Returns the .profile as
    a DataFrame (FID IID PHENO CNT CNT2 SCORESUM) or an empty frame."""
    score_file, n_dropped = dedup_score_file(Path(score_file), Path(str(out_prefix) + "_dedup.score"))
    if n_dropped:
        print(f"  ⚠️  {n_dropped} rows with duplicated chr:pos IDs (multi-allelic sites) dropped from the score file")
    r = run_plink(plink, [
        "--bfile", joint_prefix, "--score", score_file, "1", "2", "3", "sum",
        "--out", out_prefix, "--allow-extra-chr", "--threads", threads, "--memory", memory,
    ], timeout=1800)
    Path(str(out_prefix) + "_dedup.score").unlink(missing_ok=True)
    prof_path = Path(str(out_prefix) + ".profile")
    if r.returncode != 0 or not prof_path.exists():
        return pd.DataFrame()
    prof = pd.read_csv(prof_path, sep=r"\s+", dtype={"FID": str, "IID": str}, engine="python")
    score_col = "SCORESUM" if "SCORESUM" in prof.columns else "SCORE"
    prof["SCORE_VALUE"] = prof[score_col].astype(float)
    return prof


def population_distributions(ref_scores: pd.DataFrame) -> dict:
    dist = {}
    for sp in SUPER_POPS:
        s = ref_scores.loc[ref_scores["super_pop"] == sp, "SCORE_VALUE"]
        if len(s) < MIN_REF_PER_POP:
            continue
        dist[sp] = {
            "n": int(len(s)), "mean": float(s.mean()), "std": float(s.std()),
            "median": float(s.median()),
            "p5": float(np.percentile(s, 5)), "p10": float(np.percentile(s, 10)),
            "p25": float(np.percentile(s, 25)), "p75": float(np.percentile(s, 75)),
            "p90": float(np.percentile(s, 90)), "p95": float(np.percentile(s, 95)),
            "min": float(s.min()), "max": float(s.max()),
        }
    return dist


def risk_category(z: float) -> str:
    if z > 2:
        return "HIGH"
    if z > 1:
        return "ELEVATED"
    if z >= -1:
        return "AVERAGE"
    return "LOW"


# ── main ─────────────────────────────────────────────────────────────────────

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bfile", required=True, help="1000G reference PLINK prefix")
    p.add_argument("--pop-panel", required=True, help="TSV with columns sample, pop, super_pop")
    p.add_argument("--pgs-dir", required=True, help="Directory with PGS*/PGS*_clean.score files")
    p.add_argument("--user-bfile", required=True, help="User PLINK prefix (qc/qc_filtered)")
    p.add_argument("--ancestry-json", default=None, help="science/ANCESTRY_MODEL.json (assigned_population)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plink", default="plink")
    p.add_argument("--threads", default="8")
    p.add_argument("--memory", default="16000")
    p.add_argument("--min-variants-per-chrom", type=int, default=100,
                   help="A chromosome counts as covered by the user's data when it has at least this many variants")
    p.add_argument("--keep-work", action="store_true", help="Keep the joint PLINK dataset after scoring")
    p.add_argument("--site-policy", default="wgs", choices=["wgs", "array"],
                   help="wgs: sites absent from the user's data are homozygous reference (filled); "
                        "array: only sites the user was genotyped on are scored, for everyone")
    p.add_argument("--sample-prs", default=None, help=argparse.SUPPRESS)  # legacy, ignored
    args = p.parse_args(argv)

    plink = str(args.plink)
    pgs_dir = Path(args.pgs_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = out_dir / "ref_distributions"
    ref_dir.mkdir(exist_ok=True)
    if args.sample_prs:
        print("ℹ️  --sample-prs is ignored: the user is now scored jointly with the reference")

    # Population panel + reference sample IDs
    pop = pd.read_csv(args.pop_panel, sep="\t", dtype=str)
    pop_map = dict(zip(pop["sample"], pop["super_pop"]))
    print(f"Population panel: {len(pop)} samples across {pop['super_pop'].nunique()} super-pops")

    target_pop, ancestry_source = load_assigned_population(args.ancestry_json)
    print(f"Calibration population: {target_pop} ({ancestry_source})")

    # Score files + union of needed variant IDs
    pgs_files = sorted(pgs_dir.glob("PGS*/PGS*_clean.score"))
    if not pgs_files:
        sys.exit(f"❌ No PGS*/PGS*_clean.score files under {pgs_dir}")
    needed_ids = set()
    for score_file in pgs_files:
        with open(score_file) as fh:
            for line in fh:
                vid = line.split("\t", 1)[0].strip()
                if vid:
                    needed_ids.add(vid)
    print(f"{len(needed_ids):,} unique variant IDs needed across {len(pgs_files)} scores")

    # Chromosomes the user's data covers (hom-ref filling only where sequenced)
    chroms = None
    if args.site_policy == "wgs":
        chroms = user_covered_chromosomes(Path(str(args.user_bfile) + ".bim"), args.min_variants_per_chrom)
        if not chroms:
            sys.exit("❌ The user dataset covers no autosome with enough variants; cannot calibrate PGS")
        print(f"User data covers {len(chroms)} autosome(s): {', '.join(chroms)}")

    # Joint dataset
    print(f"\nBuilding joint user + reference dataset (same variants, site policy: {args.site_policy})…")
    try:
        joint, user_id_map, joint_info = build_joint_dataset(
            plink, args.bfile, args.user_bfile, needed_ids, out_dir, threads=args.threads, memory=args.memory,
            site_policy=args.site_policy, chroms=chroms, prefix="_joint")
    except JointDatasetError as e:
        sys.exit(f"❌ {e}")
    if joint_info["n_merge_conflicts_excluded"]:
        print(f"  ⚠️  {joint_info['n_merge_conflicts_excluded']} variants inconsistent between user and "
              "reference — excluded from both")
    joint_fam = read_fam_ids(Path(str(joint) + ".fam"))
    n_joint_variants = sum(1 for _ in open(str(joint) + ".bim"))
    user_iids = [iid for _, iid in joint_fam if iid in user_id_map]
    print(f"  Joint dataset: {len(joint_fam)} samples ({len(user_iids)} user), {n_joint_variants:,} variants")

    results = []
    for i, score_file in enumerate(pgs_files):
        pgs_id = score_file.parent.name
        print(f"\n[{i + 1}/{len(pgs_files)}] {pgs_id}")
        out_prefix = ref_dir / pgs_id
        prof = score_joint(plink, joint, score_file, out_prefix, args.threads, args.memory)
        if prof.empty:
            print("  ❌ No profile — skipping")
            continue

        prof["super_pop"] = prof["IID"].map(pop_map)
        ref_scores = prof[prof["super_pop"].notna() & ~prof["IID"].isin(user_id_map)]
        dist = population_distributions(ref_scores)
        with open(ref_dir / f"{pgs_id}_dist.json", "w") as f:
            json.dump(dist, f, indent=2)
        if not dist:
            print("  ❌ No population had enough reference samples — skipping")
            continue

        meta = read_score_metadata(pgs_dir, pgs_id)
        n_in_score = meta["n_snps"]
        cnt_col = "CNT" if "CNT" in prof.columns else None

        for _, row in prof[prof["IID"].isin(user_id_map)].iterrows():
            sample_score = float(row["SCORE_VALUE"])
            n_matched = int(row[cnt_col]) // 2 if cnt_col else 0
            coverage = (n_matched / n_in_score) if n_in_score else float("nan")

            pop_used, source = target_pop, ancestry_source
            if pop_used not in dist:
                print(f"  ⚠️  {pop_used} has no distribution for {pgs_id}; using EUR")
                pop_used, source = "EUR", "fallback"
            if pop_used not in dist:
                pop_used = next(iter(dist))
            d = dist[pop_used]
            z = (sample_score - d["mean"]) / d["std"] if d["std"] > 0 else 0.0
            pctl = float(scipy.stats.norm.cdf(z) * 100)
            ref_vals = ref_scores.loc[ref_scores["super_pop"] == pop_used, "SCORE_VALUE"].to_numpy()
            pctl_emp = float((ref_vals < sample_score).mean() * 100) if len(ref_vals) else float("nan")
            z_by_pop = {sp: round((sample_score - dd["mean"]) / dd["std"], 3) if dd["std"] > 0 else 0.0
                        for sp, dd in dist.items()}
            reliable = bool(n_in_score <= MAX_RELIABLE_SNPS and n_in_score > 0 and coverage >= MIN_COVERAGE)

            results.append({
                "pgs_id": pgs_id,
                "individual_id": user_id_map[row["IID"]],
                "trait": meta["trait"],
                "n_snps": n_in_score,
                "n_snps_matched": n_matched,
                "coverage": round(coverage, 4) if n_in_score else np.nan,
                "sample_score": sample_score,
                "reference_population": pop_used,
                "ancestry_source": source,
                "ref_mean": d["mean"],
                "ref_std": d["std"],
                "z_score": round(z, 3),
                "percentile": round(pctl, 1),
                "percentile_empirical": round(pctl_emp, 1),
                "risk_category": risk_category(z),
                "reliable": reliable,
                "n_populations": len(dist),
                "z_by_population": json.dumps(z_by_pop),
            })
            print(f"  {user_id_map[row['IID']]}: z={z:+.2f} vs {pop_used} "
                  f"(pctl {pctl:.1f}, matched {n_matched:,}/{n_in_score:,}, "
                  f"{'reliable' if reliable else 'unreliable'})")

        for ext in (".profile", ".log", ".nosex", ".nopred"):
            Path(str(out_prefix) + ext).unlink(missing_ok=True)

    if not args.keep_work:
        remove_joint_dataset(joint, extra_prefixes=(out_dir / "_joint_ref",))
        (out_dir / "_joint_needed.txt").unlink(missing_ok=True)

    if not results:
        sys.exit("❌ No PGS score could be calibrated")

    cal = pd.DataFrame(results).sort_values(["individual_id", "z_score"], ascending=[True, False])
    cal.to_csv(out_dir / "pgs_calibrated.csv", index=False)

    first = cal["individual_id"].iloc[0]
    cal_first = cal[cal["individual_id"] == first]
    print(f"\n{'=' * 88}")
    print(f"{'PGS ID':<12} {'Z':>8} {'Pctl':>7} {'Risk':>9} {'Matched':>18} {'Ref':>5}  Trait")
    print("-" * 88)
    for _, r in cal_first.iterrows():
        label = {"HIGH": "🔴", "ELEVATED": "🟠", "AVERAGE": "🟡"}.get(r["risk_category"], "🟢")
        print(f"{r['pgs_id']:<12} {r['z_score']:>+8.2f} {r['percentile']:>6.1f}% {label} {r['risk_category']:<8}"
              f" {r['n_snps_matched']:>8,}/{r['n_snps']:<8,} {r['reference_population']:>5}  {r['trait'][:30]}")
    n_samples = cal["individual_id"].nunique()
    suffix = "" if n_samples == 1 else f" for {n_samples} samples"
    print(f"\n✅ {len(cal_first)} PGS scores calibrated against 1000G {target_pop} ({ancestry_source}){suffix}")
    print(f"📁 Results: {out_dir / 'pgs_calibrated.csv'}")

    report = {
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        "methodology": {
            "reference_panel": f"1000 Genomes Phase 3 ({len(pop)} samples)",
            "populations": SUPER_POPS,
            "reference_population": target_pop,
            "ancestry_source": ancestry_source,
            "calibration_method": "Joint user+reference PLINK --score sum on one variant set; "
                                  "user sites absent from the VCF filled as homozygous reference; "
                                  "z-score vs the inferred super-population",
            "chromosomes_used": chroms if chroms else "all genotyped sites (array policy)",
            "site_policy": args.site_policy,
            "coverage_threshold": MIN_COVERAGE,
            "max_reliable_snps": MAX_RELIABLE_SNPS,
            "note": f"Scores with >{MAX_RELIABLE_SNPS:,} SNPs or <{MIN_COVERAGE:.0%} of their variants "
                    "in the joint set are flagged unreliable",
        },
        "individual_id": first,
        "summary": {
            "total_scores": int(len(cal_first)),
            "reliable_scores": int(cal_first["reliable"].sum()),
            "high_risk": int((cal_first["risk_category"] == "HIGH").sum()),
            "elevated_risk": int((cal_first["risk_category"] == "ELEVATED").sum()),
            "low_risk": int((cal_first["risk_category"] == "LOW").sum()),
        },
        "high_risk_traits": cal_first[cal_first["risk_category"] == "HIGH"].to_dict("records"),
        "elevated_risk_traits": cal_first[cal_first["risk_category"] == "ELEVATED"].to_dict("records"),
        "low_risk_traits": cal_first[cal_first["risk_category"] == "LOW"].to_dict("records"),
        "all_entries": cal_first.to_dict("records"),
        "other_samples": {sid: cal[cal["individual_id"] == sid].to_dict("records")
                          for sid in cal["individual_id"].unique() if sid != first},
        "reference_distributions_path": str(ref_dir) + "/",
    }
    with open(out_dir / "pgs_calibration_report.json", "w") as f:
        json.dump(report, f, indent=2, default=float)
    print(f"📁 Report: {out_dir / 'pgs_calibration_report.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
