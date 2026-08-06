#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 6 — MODULE 5: REAL PCA ADJUSTMENT                                   ║
║   scripts/pca_adjust_v2.py                                                  ║
║                                                                            ║
║   Thin CLI wrapper - the real logic lives in bluegen.ancestry              ║
║   (IMPROVEMENT_PLAN.md 2.1). This file only handles argparse/logging.       ║
║                                                                            ║
║   Output:                                                                   ║
║     prs/pca_adjusted_scores.csv                                             ║
║     prs/pca_adjustment_metrics.json                                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from bluegen.ancestry import (  # noqa: E402
    PCAdjustmentResult, PCAdjustmentReport, PCAAdjustmentV2,
)

logger = logging.getLogger(__name__)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 6 Module 5: Real PCA Adjustment (active PRS correction)"
    )
    parser.add_argument("--prs-data", required=True,
                       help="Raw PRS CSV")
    parser.add_argument("--sample-pcs", required=True,
                       help="Sample PC coordinates from pca_true_projection.py")
    parser.add_argument("--output-dir", "-o", default="prs")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--n-pcs", type=int, default=10)
    parser.add_argument("--ref-betas", help="Pre-computed reference PC betas JSON")
    parser.add_argument("--compute-ref", action="store_true",
                       help="Compute reference betas from 1000G PRS data")
    parser.add_argument("--ref-prs", help="1000G PRS data for reference betas")
    parser.add_argument("--ref-pcs", help="1000G eigenvec for reference betas")
    parser.add_argument("--population-panel", help="1000G population panel")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    adjuster = PCAAdjustmentV2(n_pcs=args.n_pcs)

    if args.compute_ref and args.ref_prs and args.ref_pcs and args.population_panel:
        adjuster.compute_ref_betas(
            ref_prs_data=args.ref_prs,
            ref_pcs_path=args.ref_pcs,
            population_panel=args.population_panel,
            output_dir=args.output_dir,
        )

    report = adjuster.adjust(
        prs_data=args.prs_data,
        sample_pcs=args.sample_pcs,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        ref_beta_path=args.ref_betas,
    )

    print(f"\n═══ PCA Adjustment Results ═══")
    print(f"  Mean |delta|: {report.mean_delta_pct:.1f}%")
    print(f"  Traits significantly adjusted: {report.traits_significantly_adjusted}/{len(report.results)}")
    print(f"\n  {'Trait':<30} {'Raw':>8} {'Adj':>8} {'Delta':>8} {'Δ%':>8}")
    print(f"  {'-'*62}")
    for r in report.results:
        sig = "*" if r.is_significant else " "
        print(f"  {r.trait[:28]:<30} {r.raw_prs:>8.3f} {r.adjusted_prs:>8.3f} "
              f"{r.delta:>+8.3f} {r.delta_pct:>+7.1f}%{sig}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
