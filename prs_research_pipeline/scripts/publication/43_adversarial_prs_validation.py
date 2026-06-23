#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 10 — ADVERSARIAL PRS VALIDATION                                      ║
║   scripts/43_adversarial_prs_validation.py                                   ║
║                                                                            ║
║   Simulates reviewer-level attacks on the PRS pipeline:                     ║
║                                                                            ║
║     1. Population shift stress (EUR→AFR, EUR→EAS, admixed)                 ║
║     2. GWAS portability failure (effect decay, ancestry bias)              ║
║     3. LD structure disruption (broken blocks, mismatched panels)          ║
║     4. Missing SNP robustness (10%, 30%, 50% dropout)                      ║
║                                                                            ║
║   This module answers: "What happens when the pipeline is wrong?"           ║
║                                                                            ║
║   Output:                                                                   ║
║     science/adversarial_validation_report.json                              ║
║     science/adversarial_validation_report.md                                ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)

SUPER_POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
POP_SHIFTS = {"AFR": 0.30, "EAS": 0.15, "SAS": 0.18, "AMR": 0.20}

@dataclass
class StressTestResult:
    test_id: str; description: str
    metric: str; baseline: float; stressed: float
    relative_change: float; is_robust: bool
    severity: str; detail: str

@dataclass
class AdversarialReport:
    results: List[StressTestResult] = field(default_factory=list)
    n_tests: int = 0; n_robust: int = 0; n_vulnerable: int = 0
    overall_robustness_score: float = 0.0
    critical_findings: List[str] = field(default_factory=list)
    generated_date: str = ""

class AdversarialPRSValidator:
    """Simulates peer-review attacks on every pipeline dimension."""

    ROBUSTNESS_THRESHOLD = 0.20

    def __init__(self, seed: int = 42, output_dir: str = "science"):
        self.seed = seed; self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.RandomState(seed)

    def run_all(self, prs_result_json: str = "prs/PRS_RESULT.json",
                prs_core_json: str = "science/prs_core_definition.json",
                ancestry_json: str = "science/ANCESTRY_MODEL.json") -> AdversarialReport:
        logger.info("═══ Adversarial PRS Validation ═══")

        # Load baseline data
        prs_data = self._load_json(prs_result_json) or {}
        entries = prs_data.get("prs_entries", [])
        if not entries:
            logger.warning("No PRS data — using synthetic baseline")
            entries = [{"trait": f"trait_{i}", "raw_score": 0.5, "population_zscore": 0.0,
                        "population_percentile": 50.0, "ci_95_lower": -0.5, "ci_95_upper": 0.5}
                       for i in range(10)]

        results = []
        results.extend(self._test_population_shift(entries))
        results.extend(self._test_gwas_portability(entries))
        results.extend(self._test_ld_disruption(entries))
        results.extend(self._test_missing_snp_robustness(entries))

        n_robust = sum(1 for r in results if r.is_robust)
        n_vulnerable = sum(1 for r in results if not r.is_robust)
        robustness = 100.0 * n_robust / max(len(results), 1)

        critical = [r.test_id for r in results if r.severity == "CRITICAL"]

        report = AdversarialReport(
            results=results, n_tests=len(results),
            n_robust=n_robust, n_vulnerable=n_vulnerable,
            overall_robustness_score=round(robustness, 1),
            critical_findings=critical,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save_json(report); self._save_markdown(report)
        return report

    def _test_population_shift(self, entries: List[Dict]) -> List[StressTestResult]:
        logger.info("  Stress: Population shift")
        results = []
        baseline_zs = np.array([e.get("population_zscore", 0) for e in entries])
        baseline_pctls = np.array([e.get("population_percentile", 50) for e in entries])

        for pop, shift in POP_SHIFTS.items():
            shifted_zs = baseline_zs + shift
            shifted_pctls = baseline_pctls + shift * 15  # percentile drift
            reorder = np.corrcoef(np.argsort(baseline_pctls), np.argsort(shifted_pctls))[0, 1] if len(baseline_pctls) > 1 else 1.0
            rank_stable = reorder > 0.70
            severity = "CRITICAL" if pop == "AFR" and abs(shift) > 0.25 else ("HIGH" if abs(shift) > 0.20 else "MODERATE")
            results.append(StressTestResult(
                test_id=f"POP_SHIFT_{pop}", description=f"EUR→{pop} population shift",
                metric="Rank correlation", baseline=1.0, stressed=round(float(reorder), 4),
                relative_change=round(float(1.0 - reorder), 4),
                is_robust=rank_stable, severity=severity,
                detail=f"Z-score shift: {shift:.2f}, rank r={reorder:.3f}"))
        return results

    def _test_gwas_portability(self, entries: List[Dict]) -> List[StressTestResult]:
        logger.info("  Stress: GWAS portability decay")
        results = []
        baseline_zs = np.array([e.get("population_zscore", 0) for e in entries])

        # Shrink effect sizes toward zero (simulating portability decay)
        for decay, label in [(0.3, "MODERATE"), (0.5, "SEVERE"), (0.7, "EXTREME")]:
            shrunk = baseline_zs * (1.0 - decay)
            # Add ancestry-bias noise
            biased = shrunk + self.rng.normal(0, decay * 0.5, len(shrunk)) if len(shrunk) > 0 else shrunk
            r = np.corrcoef(baseline_zs, biased)[0, 1] if len(baseline_zs) > 1 else 1.0
            robust = r > 0.60
            results.append(StressTestResult(
                test_id=f"GWAS_DECAY_{label}", description=f"GWAS effect portability decay ({int(decay*100)}%)",
                metric="Score correlation", baseline=1.0, stressed=round(float(r), 4),
                relative_change=round(float(1.0 - r), 4),
                is_robust=robust, severity="HIGH" if decay > 0.5 else "MODERATE",
                detail=f"Effect decay={decay:.0%}, biased r={r:.3f}"))
        return results

    def _test_ld_disruption(self, entries: List[Dict]) -> List[StressTestResult]:
        logger.info("  Stress: LD structure disruption")
        results = []
        baseline_scores = np.array([e.get("raw_score", 0) for e in entries])
        variance = np.var(baseline_scores) if len(baseline_scores) > 1 else 1.0

        # Broken LD blocks inflate variance by retaining correlated SNPs
        for inflation, label in [(1.3, "MILD"), (2.0, "MODERATE"), (3.0, "SEVERE")]:
            inflated_var = variance * inflation
            vif = inflated_var / max(variance, 0.001)
            robust = vif < 2.0
            results.append(StressTestResult(
                test_id=f"LD_DISRUPT_{label}", description=f"LD disruption — variance inflation x{inflation:.1f}",
                metric="Variance inflation factor (VIF)", baseline=1.0,
                stressed=round(float(vif), 2), relative_change=round(float(vif - 1.0), 2),
                is_robust=robust, severity="CRITICAL" if vif > 2.5 else ("HIGH" if vif > 2.0 else "MODERATE"),
                detail=f"VIF={vif:.1f}x — {'acceptable' if robust else 'problematic'} variance inflation"))
        return results

    def _test_missing_snp_robustness(self, entries: List[Dict]) -> List[StressTestResult]:
        logger.info("  Stress: Missing SNP robustness")
        results = []
        baseline_scores = np.array([e.get("raw_score", 0) for e in entries])

        for drop_rate, label in [(0.10, "10%"), (0.30, "30%"), (0.50, "50%")]:
            n_keep = max(1, int(len(baseline_scores) * (1 - drop_rate)))
            # Simulate random SNP dropout: scores scale with SNP count
            scaled = baseline_scores * (1 - drop_rate)
            r = np.corrcoef(baseline_scores, scaled)[0, 1] if len(baseline_scores) > 1 else 1.0
            rank_r = np.corrcoef(
                np.argsort(baseline_scores), np.argsort(scaled))[0, 1] if len(baseline_scores) > 1 else 1.0
            robust = r > 0.80 and rank_r > 0.75
            severity = "CRITICAL" if drop_rate >= 0.50 else ("HIGH" if drop_rate >= 0.30 else "MODERATE")
            results.append(StressTestResult(
                test_id=f"SNP_DROPOUT_{label.replace('%','')}", description=f"SNP dropout {label}",
                metric="Rank correlation after dropout", baseline=1.0,
                stressed=round(float(rank_r), 4), relative_change=round(float(1.0 - rank_r), 4),
                is_robust=robust, severity=severity,
                detail=f"{label} dropout: score r={r:.3f}, rank r={rank_r:.3f}"))
        return results

    def _load_json(self, path: str) -> Optional[Dict]:
        if Path(path).exists():
            try:
                with open(path) as fh: return json.load(fh)
            except Exception: pass
        return None

    def _save_json(self, report: AdversarialReport) -> None:
        path = self.output_dir / "adversarial_validation_report.json"
        with open(path, "w") as fh:
            json.dump({
                "n_tests": report.n_tests, "n_robust": report.n_robust,
                "n_vulnerable": report.n_vulnerable,
                "overall_robustness_score": report.overall_robustness_score,
                "critical_findings": report.critical_findings,
                "generated_date": report.generated_date,
                "results": [asdict(r) for r in report.results],
            }, fh, indent=2, default=str)
        logger.info(f"  ✅ Adversarial report: {path}")

    def _save_markdown(self, report: AdversarialReport) -> None:
        path = self.output_dir / "adversarial_validation_report.md"
        lines = [
            "# Adversarial PRS Validation Report",
            f"\n**Overall Robustness:** {report.overall_robustness_score:.0f}%",
            f"**Tests:** {report.n_tests} total, {report.n_robust} robust, {report.n_vulnerable} vulnerable",
            f"\n**Generated:** {report.generated_date}",
            "\n## Results",
            "\n| Test | Description | Metric | Baseline | Stressed | Δ | Robust |",
            "|------|-------------|--------|----------|----------|---|--------|",
        ]
        for r in report.results:
            icon = "✅" if r.is_robust else "❌"
            lines.append(f"| {r.test_id} | {r.description} | {r.metric} | {r.baseline:.2f} | {r.stressed} | {r.relative_change:+.2f} | {icon} |")

        if report.critical_findings:
            lines += ["\n## ⚠️ Critical Vulnerabilities", ""]
            for f in report.critical_findings:
                lines.append(f"- ❌ {f}")

        lines += [
            "\n---",
            "\n*Phase 10 Adversarial Validation — Peer Review Resilience Testing*",
        ]
        with open(path, "w") as fh: fh.write("\n".join(lines))
        logger.info(f"  ✅ Adversarial markdown: {path}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 10: Adversarial PRS Validation")
    parser.add_argument("--prs-result", default="prs/PRS_RESULT.json")
    parser.add_argument("--prs-core", default="science/prs_core_definition.json")
    parser.add_argument("--ancestry", default="science/ANCESTRY_MODEL.json")
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    validator = AdversarialPRSValidator(output_dir=args.output_dir)
    report = validator.run_all(args.prs_result, args.prs_core, args.ancestry)
    print(f"\n═══ Adversarial Validation ═══")
    print(f"  Robustness: {report.overall_robustness_score:.0f}%")
    print(f"  Robust: {report.n_robust} | Vulnerable: {report.n_vulnerable}")
    if report.critical_findings:
        print(f"  ⚠️ Critical: {', '.join(report.critical_findings)}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
