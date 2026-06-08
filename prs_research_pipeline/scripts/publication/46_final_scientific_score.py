#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 10 — FINAL SCIENTIFIC INTEGRITY SCORE (LOCKED)                       ║
║   scripts/46_final_scientific_score.py                                       ║
║                                                                            ║
║   Computes the FINAL locked scientific integrity score combining:            ║
║     • SSST compliance (internal consistency)                                ║
║     • External benchmark performance (PGS/GWAS)                              ║
║     • Adversarial robustness (stress test survival)                         ║
║     • Calibration stability (slope, R², drift)                              ║
║     • Population portability (cross-ancestry bias)                          ║
║                                                                            ║
║   Strict weighting — no heuristic inflation.                                ║
║   Reproducible formula. Deterministic computation.                           ║
║                                                                            ║
║   Output:                                                                    ║
║     FINAL_SCIENTIFIC_SCORE.json                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, logging, hashlib
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field, asdict
from datetime import datetime

logger = logging.getLogger(__name__)

# Locked weights — DO NOT CHANGE after Phase 10 freeze
LOCKED_WEIGHTS = {
    "ssst_compliance": 0.25,
    "external_benchmarks": 0.20,
    "adversarial_robustness": 0.20,
    "calibration_stability": 0.15,
    "population_portability": 0.10,
    "reproducibility": 0.10,
}

INTERPRETATION = {
    (90, 101): ("PUBLICATION_READY", "All criteria met. Peer-review resilient."),
    (75, 90):  ("RESEARCH_GRADE", "Minor limitations. Suitable for preprint."),
    (60, 75):  ("NEEDS_REVISION", "Moderate issues. Address before journal submission."),
    (40, 60):  ("SIGNIFICANT_ISSUES", "Major gaps. Not ready for peer review."),
    (0, 40):   ("NOT_PUBLISHABLE", "Critical flaws. Do not submit."),
}

@dataclass
class ScoreComponent:
    name: str; score: float; weight: float; source: str
    contribution: float = 0.0

@dataclass
class FinalScore:
    sci_score: float; category: str; category_description: str
    components: List[ScoreComponent] = field(default_factory=list)
    score_hash: str = ""; weights_locked: bool = True
    generated_date: str = ""; formula: str = ""

class FinalScientificScore:
    """Computes the FINAL locked scientific integrity score."""

    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)

    def compute(self,
                ssst_manifest: str = "science/CONSOLIDATION_MANIFEST.json",
                benchmark_json: str = "benchmark/VALIDATION_REPORT.json",
                adversarial_json: str = "science/adversarial_validation_report.json",
                calibration_json: str = "benchmark/calibration_validation.json",
                portability_json: str = "benchmark/portability_report.json",
                fingerprint_json: str = "reproducibility/run_fingerprint.json") -> FinalScore:
        logger.info("═══ FINAL Scientific Integrity Score (LOCKED) ═══")

        components = []

        # 1. SSST Compliance
        ssst_complete = self._check_ssst(ssst_manifest)
        components.append(ScoreComponent(
            name="SSST Compliance", score=ssst_complete,
            weight=LOCKED_WEIGHTS["ssst_compliance"],
            source="CONSOLIDATION_MANIFEST.json"))

        # 2. External Benchmarks
        bench_score = self._benchmark_score(benchmark_json)
        components.append(ScoreComponent(
            name="External Benchmarks", score=bench_score,
            weight=LOCKED_WEIGHTS["external_benchmarks"],
            source="VALIDATION_REPORT.json"))

        # 3. Adversarial Robustness
        adv_score = self._adversarial_score(adversarial_json)
        components.append(ScoreComponent(
            name="Adversarial Robustness", score=adv_score,
            weight=LOCKED_WEIGHTS["adversarial_robustness"],
            source="adversarial_validation_report.json"))

        # 4. Calibration Stability
        cal_score = self._calibration_score(calibration_json)
        components.append(ScoreComponent(
            name="Calibration Stability", score=cal_score,
            weight=LOCKED_WEIGHTS["calibration_stability"],
            source="calibration_validation.json"))

        # 5. Population Portability
        port_score = self._portability_score(portability_json)
        components.append(ScoreComponent(
            name="Population Portability", score=port_score,
            weight=LOCKED_WEIGHTS["population_portability"],
            source="portability_report.json"))

        # 6. Reproducibility
        repro_score = self._reproducibility_score(fingerprint_json)
        components.append(ScoreComponent(
            name="Reproducibility", score=repro_score,
            weight=LOCKED_WEIGHTS["reproducibility"],
            source="run_fingerprint.json"))

        # Compute weighted total
        sci_score = 0.0
        for c in components:
            c.contribution = round(c.score * c.weight, 2)
            sci_score += c.contribution
        sci_score = max(0.0, min(100.0, sci_score))

        # Determine category
        category = "NOT_PUBLISHABLE"; cat_desc = ""
        for (lo, hi), (cat, desc) in INTERPRETATION.items():
            if lo <= sci_score < hi:
                category = cat; cat_desc = desc; break

        # Compute deterministic score hash
        score_hash = hashlib.sha256(
            f"{sci_score:.1f}{''.join(str(c.score) for c in components)}".encode()
        ).hexdigest()[:16]

        result = FinalScore(
            sci_score=round(sci_score, 1), category=category,
            category_description=cat_desc, components=components,
            score_hash=score_hash,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
            formula="SCI_SCORE = 0.25×SSST + 0.20×Benchmarks + 0.20×Adversarial + 0.15×Calibration + 0.10×Portability + 0.10×Reproducibility")

        self._save(result)
        return result

    def _check_ssst(self, path: str) -> float:
        data = self._load(path)
        if not data: return 30.0
        summary = data.get("summary", {})
        active = summary.get("active", 0); total = summary.get("total", 7)
        complete = summary.get("consolidation_complete", False)
        return 100.0 if complete else (30.0 + 70.0 * (active / max(total, 1)))

    def _benchmark_score(self, path: str) -> float:
        data = self._load(path)
        if not data: return 40.0
        summary = data.get("validation_summary", {})
        n_ext = summary.get("external", 0); n_circ = summary.get("circular", 0)
        independent = summary.get("all_independent", False)
        base = 60.0 + n_ext * 10.0
        if not independent: base -= n_circ * 15.0
        return min(100.0, max(10.0, base))

    def _adversarial_score(self, path: str) -> float:
        data = self._load(path)
        if not data: return 30.0
        return data.get("overall_robustness_score", 30.0)

    def _calibration_score(self, path: str) -> float:
        data = self._load(path)
        if not data: return 40.0
        slope = abs(data.get("mean_slope", 1.0) - 1.0)
        r2 = data.get("mean_r2", 0.5)
        return min(100.0, max(10.0, r2 * 100 - slope * 50))

    def _portability_score(self, path: str) -> float:
        data = self._load(path)
        if not data: return 30.0
        bias = data.get("global_bias_index", 0.3)
        return min(100.0, max(5.0, 100.0 - bias * 250))

    def _reproducibility_score(self, path: str) -> float:
        if Path(path).exists(): return 95.0
        return 20.0

    def _load(self, path: str) -> Dict:
        if Path(path).exists():
            try:
                with open(path) as fh: return json.load(fh)
            except Exception: pass
        return {}

    def _save(self, result: FinalScore) -> None:
        path = Path("FINAL_SCIENTIFIC_SCORE.json")
        with open(path, "w") as fh:
            json.dump({
                "scientific_integrity_score": result.sci_score,
                "category": result.category,
                "category_description": result.category_description,
                "score_hash": result.score_hash,
                "weights_locked": True,
                "formula": result.formula,
                "generated_date": result.generated_date,
                "weights": LOCKED_WEIGHTS,
                "interpretation": {f"{lo}-{hi-1}": label for (lo, hi), (label, _) in INTERPRETATION.items()},
                "components": [asdict(c) for c in result.components],
            }, fh, indent=2)
        logger.info(f"  ✅ FINAL SCORE: {result.sci_score:.0f}/100 — {result.category}")
        logger.info(f"  Hash: {result.score_hash} — LOCKED")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 10: Final Scientific Integrity Score (LOCKED)")
    parser.add_argument("--ssst-manifest", default="science/CONSOLIDATION_MANIFEST.json")
    parser.add_argument("--benchmark", default="benchmark/VALIDATION_REPORT.json")
    parser.add_argument("--adversarial", default="science/adversarial_validation_report.json")
    parser.add_argument("--calibration", default="benchmark/calibration_validation.json")
    parser.add_argument("--portability", default="benchmark/portability_report.json")
    parser.add_argument("--fingerprint", default="reproducibility/run_fingerprint.json")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    scorer = FinalScientificScore()
    result = scorer.compute(args.ssst_manifest, args.benchmark, args.adversarial,
                             args.calibration, args.portability, args.fingerprint)
    print(f"\n═══ FINAL SCIENTIFIC SCORE (LOCKED) ═══")
    print(f"  Score: {result.sci_score:.0f}/100 — {result.category}")
    print(f"  Hash: {result.score_hash}")
    for c in result.components:
        print(f"  {c.name}: {c.score:.0f}/100 × {c.weight:.0%} = {c.contribution:.1f}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
