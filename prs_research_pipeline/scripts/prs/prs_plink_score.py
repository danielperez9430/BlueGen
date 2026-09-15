#!/usr/bin/env python3
"""
PRS Computation via PLINK --score — multi-trait, multi-sample.

Thin CLI wrapper - the real logic lives in bluegen.scoring
(IMPROVEMENT_PLAN.md 2.1). This file only handles argparse.

Computes polygenic risk scores using the weighted sum method:
    PRS_i = Σ (β_j × G_ij)

Input:
    --snp-db       CSV with columns: chrom, pos, effect_allele, weight, trait_category
    --bfile        PLINK binary dataset prefix (qc/qc_filtered)
    --output-dir   Output directory (default: prs/)
    --plink        Path to PLINK binary (default: plink)

Output:
    prs/prs_raw.csv   Columns: individual_id, trait, prs_raw, n_snps, n_snps_used

Usage:
    python prs_plink_score.py \\
        --snp-db data/snp_database_annotated.csv \\
        --bfile qc/qc_filtered \\
        --output-dir prs/ \\
        --plink tools/plink
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from bluegen.scoring import (  # noqa: E402
    fix_duplicate_bim_ids, build_score_rows, write_score_file,
    run_plink_score, parse_plink_profile, compute_prs_plink_score, compute_prs_joint,
)


def main():
    parser = argparse.ArgumentParser(description="PRS Computation via PLINK --score")
    parser.add_argument("--snp-db", required=True, help="Path to SNP database CSV")
    parser.add_argument("--bfile", required=True, help="PLINK bfile prefix (qc/qc_filtered)")
    parser.add_argument("--output-dir", default="prs/", help="Output directory")
    parser.add_argument("--plink", default="plink", help="Path to PLINK binary")
    parser.add_argument("--threads", type=int, default=4, help="PLINK threads")
    parser.add_argument("--memory", type=int, default=8000, help="PLINK memory (MB)")
    # Joint mode (RELEASE_PLAN 3.0.3): score the user together with the 1000G
    # reference on one variant set; also writes prs_reference_raw.csv.
    parser.add_argument("--ref-bfile", help="1000G PLINK prefix → joint scoring (recommended)")
    parser.add_argument("--pop-panel", help="1000G population panel (sample pop super_pop)")
    parser.add_argument("--site-policy", default="wgs", choices=["wgs", "array"],
                        help="wgs: absent sites are homozygous reference; array: only genotyped sites")
    args = parser.parse_args()

    if args.ref_bfile:
        compute_prs_joint(
            snp_db=args.snp_db,
            user_bfile=args.bfile,
            ref_bfile=args.ref_bfile,
            pop_panel=args.pop_panel,
            output_dir=args.output_dir,
            plink=args.plink,
            site_policy=args.site_policy,
            threads=args.threads,
            memory=args.memory,
        )
    else:
        print("  PRS: user-only scoring (no --ref-bfile) — absent panel SNPs are dropped, "
              "not treated as homozygous reference; z-scores vs precomputed distributions are biased")
        compute_prs_plink_score(
            snp_db=args.snp_db,
            bfile=args.bfile,
            output_dir=args.output_dir,
            plink=args.plink,
            threads=args.threads,
            memory=args.memory,
        )


if __name__ == "__main__":
    main()
