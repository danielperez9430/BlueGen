#!/usr/bin/env python3
"""
Build population reference distributions from 1000 Genomes + curated SNP database.
One-time computation — cached forever after.

Runs PLINK --score on genome-wide 1000G reference, aggregates PRS
per trait × population, and saves empirical distribution parameters.
"""

import sys
import os
import json
import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SUPER_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

def build_reference_distributions(
    plink_bin: str = "plink",
    ref_bfile: str = "reference/1000G_full/1000G_full",
    snp_db: str = "data/snp_database_annotated.csv",
    pop_panel: str = "reference/1000G_full/population_panel.txt",
    output_dir: str = "reference/population_distributions",
    threads: int = 4,
    memory: int = 16000,
):
    logger.info("═══ Building Population Reference Distributions ═══")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load SNP database
    db = pd.read_csv(snp_db, dtype=str)
    trait_col = "trait_category" if "trait_category" in db.columns else "trait"
    traits = sorted(db[trait_col].dropna().unique())
    logger.info(f"  Traits: {len(traits)}")

    # Load population panel
    panel = pd.read_csv(pop_panel, sep=r"\s+", dtype=str)
    sample_to_pop = {}
    for _, row in panel.iterrows():
        sample_to_pop[str(row.iloc[0])] = str(row.iloc[2])

    all_scores = []
    hwe_fallbacks = []

    for trait in traits:
        trait_snps = db[db[trait_col] == trait]
        if len(trait_snps) == 0:
            continue

        # Create PLINK score file for this trait
        # Format: variant_id allele weight
        # Build variant_id from chrom:pos to match BIM format (no "chr" prefix)
        score_lines = []
        seen_ids = set()
        for _, row in trait_snps.iterrows():
            chrom = str(row.get("chrom", "")).replace("chr", "")
            pos = row.get("pos", "")
            allele = row.get("effect_allele", "")
            weight = row.get("weight", row.get("beta", "1.0"))
            direction = str(row.get("effect_direction", "+")).strip()
            if chrom and pos and allele:
                try:
                    w = float(weight)
                    if direction == "-":
                        w = -w
                    var_id = f"{chrom}:{pos}"
                    # Deduplicate: same position, different rsID (multi-allelic/co-located)
                    if var_id in seen_ids:
                        continue
                    seen_ids.add(var_id)
                    score_lines.append(f"{var_id}\t{allele}\t{w}")
                except ValueError:
                    continue

        if not score_lines:
            logger.warning(f"  {trait}: no valid SNPs for scoring")
            continue

        # Write temporary score file (sanitize trait name for filesystem)
        safe_trait = trait.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '').replace('&', 'and')
        score_file = output_dir / f"_{safe_trait}.score"
        with open(score_file, "w") as fh:
            fh.write("\n".join(score_lines))

        # Run PLINK --score on 1000G reference
        profile_file = output_dir / f"_{trait.replace(' ', '_')}.profile"
        cmd = [
            plink_bin, "--bfile", ref_bfile,
            "--score", str(score_file), "1", "2", "3",
            "--out", str(output_dir / f"_{trait.replace(' ', '_')}"),
            "--threads", str(threads), "--memory", str(memory),
            "--allow-extra-chr",
        ]
        logger.info(f"  Scoring {trait} ({len(score_lines)} SNPs)...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            logger.warning(f"    PLINK score failed for {trait} — HWE fallback")
            # Fallback: Hardy-Weinberg population estimates
            fallback_mu = 0.0
            fallback_var = 0.0
            for line in score_lines:
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        beta = float(parts[2])
                        p = 0.25
                        fallback_mu += 2 * p * beta
                        fallback_var += 2 * p * (1-p) * beta * beta
                    except ValueError:
                        continue
            sigma = (fallback_var ** 0.5) if fallback_var > 0 else 0.01
            # Store for later merging into distributions dict
            hwe_fallbacks.append((trait, fallback_mu, sigma))
            logger.info(f"    μ={fallback_mu:.4f}, σ={sigma:.4f}")
            score_file.unlink(missing_ok=True)
            continue

        # Read profile and aggregate
        if profile_file.exists():
            try:
                prof = pd.read_csv(profile_file, sep=r"\s+", dtype={"IID": str})
                for _, row in prof.iterrows():
                    all_scores.append({
                        "individual_id": str(row["IID"]),
                        "trait": trait,
                        "prs_raw": float(row.get("SCORE", row.get("SCORESUM", 0))),
                    })
            except Exception as e:
                logger.warning(f"    Could not parse profile for {trait}: {e}")

        # Clean temp files
        score_file.unlink(missing_ok=True)
        profile_file.unlink(missing_ok=True)
        # Also clean PLINK log/nosex
        for ext in [".log", ".nosex", ".profile"]:
            f = output_dir / f"_{trait.replace(' ', '_')}{ext}"
            f.unlink(missing_ok=True)

    if not all_scores:
        logger.error("  No scores computed — check PLINK and score files")
        return None

    # Build DataFrame
    scores_df = pd.DataFrame(all_scores)
    scores_df["population"] = scores_df["individual_id"].map(sample_to_pop)
    scores_df = scores_df[scores_df["population"].isin(SUPER_POPULATIONS)]
    n_scored = scores_df["individual_id"].nunique()
    logger.info(f"  Scored {n_scored} individuals across {len(all_scores)} trait-samples")

    # Compute per-population distributions
    distributions = {}
    for trait in sorted(scores_df["trait"].unique()):
        trait_data = scores_df[scores_df["trait"] == trait]
        distributions[trait] = {}
        for pop in SUPER_POPULATIONS:
            pop_data = trait_data[trait_data["population"] == pop]["prs_raw"]
            n = len(pop_data)
            if n < 5:
                continue
            values = pop_data.values.astype(np.float64)
            mean = float(np.mean(values))
            std = float(np.std(values, ddof=1)) if n > 1 else 1.0
            median = float(np.median(values))
            q25, q75 = np.percentile(values, [25, 75]) if n >= 4 else (0, 0)
            q25_val = float(np.percentile(values, 25)) if n >= 4 else mean - 0.674 * std
            q75_val = float(np.percentile(values, 75)) if n >= 4 else mean + 0.674 * std
            p5_val = float(np.percentile(values, 5)) if n >= 20 else mean - 1.645 * std
            p95_val = float(np.percentile(values, 95)) if n >= 20 else mean + 1.645 * std
            # Compute skewness and kurtosis if scipy available
            try:
                from scipy import stats as scipy_stats
                skew_val = float(scipy_stats.skew(values))
                kurt_val = float(scipy_stats.kurtosis(values))
                _, shapiro_p = scipy_stats.shapiro(values[:min(n, 5000)])
                shapiro_p = float(shapiro_p)
            except ImportError:
                skew_val = 0.0; kurt_val = 0.0; shapiro_p = 1.0

            distributions[trait][pop] = {
                "trait": trait, "population": pop, "n_samples": int(n),
                "mean": round(mean, 6), "std": round(std if std > 0 else 1.0, 6),
                "median": round(median, 6), "iqr": round(float(q75_val - q25_val), 6),
                "percentile_5": round(p5_val, 6),
                "percentile_25": round(q25_val, 6),
                "percentile_75": round(q75_val, 6),
                "percentile_95": round(p95_val, 6),
                "skewness": round(skew_val, 6),
                "kurtosis": round(kurt_val, 6),
                "shapiro_p": round(shapiro_p, 6),
            }

    # Merge HWE fallbacks into distributions
    for trait, mu, sigma in hwe_fallbacks:
        if trait not in distributions:
            distributions[trait] = {}
        for pop in SUPER_POPULATIONS:
            distributions[trait][pop] = {
                "trait": trait, "population": pop, "n_samples": 0,
                "mean": round(mu, 6), "std": round(sigma, 6),
                "median": round(mu, 6), "q25": round(mu - 0.674*sigma, 6),
                "q75": round(mu + 0.674*sigma, 6),
                "p5": round(mu - 1.645*sigma, 6), "p95": round(mu + 1.645*sigma, 6),
                "skewness": 0.0, "kurtosis": 0.0, "shapiro_p": 1.0,
                "method": "Hardy-Weinberg estimate (SNP not in 1000G)",
            }
    logger.info(f"  HWE fallbacks merged: {len(hwe_fallbacks)} traits")

    # Save distributions
    dist_path = output_dir / "reference_distributions.json"
    with open(dist_path, "w") as fh:
        json.dump({
            "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
            "reference_panel": "1000 Genomes Phase 3 (genome-wide)",
            "n_populations": len(SUPER_POPULATIONS),
            "n_traits": len(distributions),
            "n_reference_samples": n_scored,
            "calibration_method": "Empirical 1000G population-stratified distributions",
            "distributions": distributions,
        }, fh, indent=2)

    # Also save CSV for compatibility
    csv_rows = []
    for trait, pops in distributions.items():
        for pop, d in pops.items():
            csv_rows.append(d)
    pd.DataFrame(csv_rows).to_csv(output_dir / "reference_distributions.csv", index=False)

    logger.info(f"  ✅ Distributions saved: {dist_path}")
    logger.info(f"  Traits: {len(distributions)} | Populations: {len(SUPER_POPULATIONS)}")
    return distributions


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build population reference distributions")
    parser.add_argument("--plink", default="plink", help="PLINK binary path (auto-detected if not specified)")
    parser.add_argument("--bfile", default="reference/1000G_full/1000G_full")
    parser.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    parser.add_argument("--pop-panel", default="reference/1000G_full/population_panel.txt")
    parser.add_argument("--output-dir", "-o", default="reference/population_distributions")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory", type=int, default=16000)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")

    # Auto-detect PLINK if not explicitly provided
    if args.plink == "plink":
        candidates = [
            str(Path(__file__).parent.parent.parent / "tools" / "plink"),
            str(Path(__file__).parent.parent / "tools" / "plink"),
            str(Path.cwd() / "tools" / "plink"),
        ]
        for c in candidates:
            if os.path.exists(c) and os.access(c, os.X_OK):
                args.plink = c
                logger.info(f"Auto-detected PLINK: {c}")
                break
        else:
            system_plink = shutil.which("plink") or shutil.which("plink2")
            if system_plink:
                args.plink = system_plink
            else:
                logger.error("PLINK not found in tools/ or PATH")
                return 1

    build_reference_distributions(
        plink_bin=args.plink, ref_bfile=args.bfile,
        snp_db=args.snp_db, pop_panel=args.pop_panel,
        output_dir=args.output_dir, threads=args.threads, memory=args.memory,
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
