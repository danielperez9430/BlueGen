#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 6 — MODULE 6: REAL POPULATION CALIBRATION                           ║
║   scripts/population_calibrate_v2.py                                        ║
║                                                                            ║
║   Thin CLI wrapper - the real logic lives in bluegen.calibration            ║
║   (IMPROVEMENT_PLAN.md 2.1). This file only handles argparse/logging.       ║
║                                                                            ║
║   Output:                                                                   ║
║     reference/population_distributions/                                    ║
║     prs/population_calibrated_v2.csv                                        ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from bluegen.calibration import (  # noqa: E402
    SUPER_POPULATIONS, PopulationDistribution, CalibratedPRS, PopulationCalibrationV2,
)

logger = logging.getLogger(__name__)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 6 Module 6: Real Population Calibration (empirical distributions)"
    )
    parser.add_argument("--ref-prs", help="1000G PRS data CSV (for building distributions)")
    parser.add_argument("--population-panel", help="1000G population panel")
    parser.add_argument("--sample-prs", help="Sample PRS CSV (for calibration)")
    parser.add_argument("--ancestry-json", help="Ancestry classification JSON")
    parser.add_argument("--ref-dist-dir", help="Pre-computed reference distribution directory")
    parser.add_argument("--output-dir", "-o", default="prs")
    parser.add_argument("--dist-output-dir", default="reference/population_distributions")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--build-only", action="store_true",
                       help="Only build reference distributions")
    parser.add_argument("--calibrate-only", action="store_true",
                       help="Only calibrate sample")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    calibrator = PopulationCalibrationV2()

    if args.ref_prs and args.population_panel and not args.calibrate_only:
        calibrator.build_reference_distributions(
            ref_prs=args.ref_prs,
            population_panel=args.population_panel,
            output_dir=args.dist_output_dir,
        )

    if args.sample_prs and args.ancestry_json and not args.build_only:
        calibrated = calibrator.calibrate_sample(
            sample_prs=args.sample_prs,
            ancestry_json=args.ancestry_json,
            output_dir=args.output_dir,
            sample_id=args.sample_id,
            ref_dist_dir=args.ref_dist_dir or args.dist_output_dir,
        )

        print(f"\n═══ Population-Calibrated PRS (V2) ═══")
        print(f"  {'Trait':<30} {'PRS':>8} {'z_pop':>8} {'%_pop':>8} {'Risk':>8}")
        print(f"  {'-'*62}")
        for c in calibrated:
            print(f"  {c.trait[:28]:<30} {c.prs_raw:>8.3f} {c.z_score_population:>+8.2f} "
                  f"{c.percentile_population:>7.1f}% {c.risk_category:>8}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
