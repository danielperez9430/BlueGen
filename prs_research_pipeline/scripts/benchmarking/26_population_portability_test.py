#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 8 — MODULE 4: POPULATION PORTABILITY TEST                           ║
║   scripts/26_population_portability_test.py                                 ║
║                                                                            ║
║   Evaluates PRS behavior across all five 1000 Genomes super-populations.    ║
║                                                                            ║
║   Metrics per population:                                                   ║
║     • PRS score shift (mean difference from EUR reference)                  ║
║     • Calibration drift (expected vs observed z-score distribution)         ║
║     • Rank instability (percentile re-ranking across populations)           ║
║     • Ancestry bias index (composite score of portability)                  ║
║                                                                            ║
║   Output:                                                                   ║
║     benchmark/portability_report.json                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np
import pandas as pd
try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None

logger = logging.getLogger(__name__)

SUPER_POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

@dataclass
class PortabilityMetrics:
    population: str
    mean_prs_shift: float = 0.0
    calibration_drift: float = 0.0
    rank_instability: float = 0.0
    ancestry_bias_index: float = 0.0
    n_reference_samples: int = 0
    status: str = "UNKNOWN"

@dataclass
class PortabilityReport:
    metrics: List[PortabilityMetrics] = field(default_factory=list)
    global_bias_index: float = 0.0
    most_biased_population: str = ""
    least_biased_population: str = ""
    generated_date: str = ""

class PopulationPortabilityTester:
    """Tests PRS portability across super-populations."""

    def __init__(self, output_dir: str = "benchmark"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def test(self, calibration_csv: Optional[str] = None,
             population_panel: Optional[str] = None,
             ancestry_json: Optional[str] = None,
             reference_distributions_path: Optional[str] = None) -> PortabilityReport:
        logger.info("═══ Population Portability Test ═══")

        # ── Load reference distributions for dynamic computation ────────────
        if reference_distributions_path is None:
            reference_distributions_path = "reference/population_distributions/reference_distributions.json"

        ref_dists = {}
        ref_path = Path(reference_distributions_path)
        if not ref_path.exists():
            logger.error(f"  Reference distributions not found: {reference_distributions_path}")
            raise FileNotFoundError(f"Reference distributions required: {reference_distributions_path}")

        with open(ref_path) as fh:
            ref_dists = json.load(fh)
        distributions = ref_dists.get("distributions", {})
        n_traits = len(distributions)
        n_ref_samples = ref_dists.get("n_reference_samples", 2504)
        logger.info(f"  Loaded {n_traits} traits from reference distributions ({n_ref_samples} samples)")

        # ── Compute per-population metrics dynamically ──────────────────────
        pop_shifts = {pop: [] for pop in SUPER_POPS}
        pop_drifts = {pop: [] for pop in SUPER_POPS}
        pop_ranks = {pop: [] for pop in SUPER_POPS}
        eur_ranks = []

        for trait, pops in distributions.items():
            eur = pops.get("EUR", {})
            eur_mean = eur.get("mean", 0)
            eur_std = eur.get("std", 1.0)
            eur_median = eur.get("median", eur.get("p50", 0))

            # Skip traits with default/fake distribution data
            if eur_std <= 0 or (eur_std == 1.0 and eur_mean == 0.0):
                continue

            eur_ranks.append(eur_median)

            for pop in SUPER_POPS:
                pdata = pops.get(pop, {})
                pop_mean = pdata.get("mean", 0)
                pop_std = pdata.get("std", 1.0)
                pop_median = pdata.get("median", pdata.get("p50", 0))
                if pop_std <= 0:
                    pop_std = 1.0

                # shift: Cohen's d-style (pooled SD for stability)
                pooled_sd = np.sqrt((eur_std**2 + pop_std**2) / 2)
                if pooled_sd <= 0.001:
                    pooled_sd = 0.001
                shift = min(abs(pop_mean - eur_mean) / pooled_sd, 2.0)
                pop_shifts[pop].append(shift)

                # drift: relative SD difference, capped
                drift = min(abs(pop_std - eur_std) / max(eur_std, 0.001), 2.0)
                pop_drifts[pop].append(drift)

                # rank for Spearman correlation
                pop_ranks[pop].append(pop_median)

        # Compute Spearman rank instability
        # Use scipy if available, otherwise fall back to Pearson on ranks
        rank_instability = {}
        try:
            from scipy import stats as scipy_stats
            for pop in SUPER_POPS:
                if pop == "EUR":
                    rank_instability[pop] = 0.0
                else:
                    rho, _ = scipy_stats.spearmanr(eur_ranks, pop_ranks[pop])
                    rank_instability[pop] = max(0.0, 1.0 - rho)
        except ImportError:
            # Fallback: fraction of traits where ranking differs substantially from EUR
            logger.info("  scipy not available — using rank reversal fraction for instability")
            for pop in SUPER_POPS:
                if pop == "EUR":
                    rank_instability[pop] = 0.0
                else:
                    reversals = 0
                    for pr, er in zip(pop_ranks[pop], eur_ranks):
                        if abs(er) > 0.001 and abs(pr - er) / abs(er) > 0.5:
                            reversals += 1
                    rank_instability[pop] = min(1.0, reversals / max(len(eur_ranks), 1))

        # Reference sample counts from 1000 Genomes Phase 3
        ref_counts = {"EUR": 503, "AFR": 661, "EAS": 504, "SAS": 489, "AMR": 347}

        metrics = []
        for pop in SUPER_POPS:
            # Use mean across traits with per-trait capping at 2.0 (already applied above)
            shift = min(float(np.mean(pop_shifts[pop])), 1.0) if pop_shifts[pop] else 0.0
            drift = min(float(np.mean(pop_drifts[pop])), 1.0) if pop_drifts[pop] else 0.0
            instability = min(rank_instability.get(pop, 0.0), 1.0)
            abi = (shift + drift + instability) / 3

            if abi < 0.10:
                status = "GOOD_PORTABILITY"
            elif abi < 0.25:
                status = "MODERATE_PORTABILITY"
            else:
                status = "LIMITED_PORTABILITY"

            metrics.append(PortabilityMetrics(
                population=pop,
                mean_prs_shift=round(shift, 4),
                calibration_drift=round(drift, 4),
                rank_instability=round(instability, 4),
                ancestry_bias_index=round(abi, 4),
                n_reference_samples=ref_counts.get(pop, 0),
                status=status))

        global_bias = float(np.mean([m.ancestry_bias_index for m in metrics]))
        most = max(metrics, key=lambda m: m.ancestry_bias_index)
        least = min(metrics, key=lambda m: m.ancestry_bias_index)

        report = PortabilityReport(
            metrics=metrics,
            global_bias_index=round(float(global_bias), 4),
            most_biased_population=most.population,
            least_biased_population=least.population,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save_report(report)
        return report

    def _save_report(self, report: PortabilityReport) -> None:
        path = self.output_dir / "portability_report.json"
        with open(path, "w") as fh:
            json.dump({
                "global_bias_index": report.global_bias_index,
                "most_biased": report.most_biased_population,
                "least_biased": report.least_biased_population,
                "generated_date": report.generated_date,
                "populations": [asdict(m) for m in report.metrics],
            }, fh, indent=2)
        logger.info(f"  ✅ Portability: {path}")
        logger.info(f"  Global bias index: {report.global_bias_index:.3f}")
        logger.info(f"  Most biased: {report.most_biased_population}")
        logger.info(f"  Least biased: {report.least_biased_population}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 8 M4: Population Portability Test")
    parser.add_argument("--calibration-csv", help="Population calibration CSV")
    parser.add_argument("--population-panel", help="1000G population panel")
    parser.add_argument("--ancestry-json", help="Ancestry classification")
    parser.add_argument("--reference-distributions", default="reference/population_distributions/reference_distributions.json",
                       help="Reference distributions JSON (56 traits x 5 populations)")
    parser.add_argument("--output-dir", "-o", default="benchmark")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    tester = PopulationPortabilityTester(args.output_dir)
    report = tester.test(args.calibration_csv, args.population_panel, args.ancestry_json,
                         reference_distributions_path=args.reference_distributions)
    print(f"\n═══ Population Portability ═══")
    print(f"  Global bias: {report.global_bias_index:.3f}")
    for m in report.metrics:
        print(f"  {m.population}: bias={m.ancestry_bias_index:.3f} ({m.status})")
    return 0

if __name__ == "__main__":
    sys.exit(main())
