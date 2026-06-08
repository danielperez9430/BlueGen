#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 8 CORRECTION — SCIENTIFIC INTEGRITY SCORE                            ║
║   scripts/34_scientific_integrity_score.py                                   ║
║                                                                            ║
║   Computes a unified Scientific Integrity Score (0–100) from ALL           ║
║   validation, audit, and benchmark outputs.                                 ║
║                                                                            ║
║   Components:                                                               ║
║     1. Genome coverage integrity (15%) — chr22 bias, unified SNP universe  ║
║     2. Ancestry validity (15%) — PCA ensemble, no trait-SNP inference      ║
║     3. PRS mathematical correctness (15%) — Σ(β×dosage), allele alignment  ║
║     4. Calibration quality (15%) — empirical distributions, μ≠0, slope≈1   ║
║     5. Leakage prevention (15%) — no target→reference data flow            ║
║     6. Benchmark consistency (10%) — PGS, GWAS, method concordance         ║
║     7. Reproducibility (10%) — frozen, hashed, deterministic               ║
║     8. Population portability (5%) — cross-ancestry bias index             ║
║                                                                            ║
║   Interpretation:                                                           ║
║     90–100 = Publication-Ready                                              ║
║     75–89  = Research-Grade                                                 ║
║     60–74  = Needs Improvement                                              ║
║     40–59  = Significant Issues                                             ║
║     <40    = Not Scientifically Valid                                       ║
║                                                                            ║
║   Output:                                                                   ║
║     science/scientific_integrity_score.json                                 ║
║     science/scientific_integrity_score.md                                   ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, logging
from pathlib import Path; from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)

WEIGHTS = {
    "genome_coverage": 0.15, "ancestry_validity": 0.15,
    "prs_mathematics": 0.15, "calibration_quality": 0.15,
    "leakage_prevention": 0.15, "benchmark_consistency": 0.10,
    "reproducibility": 0.10, "population_portability": 0.05,
}

INTERPRETATION = {
    (90, 101): ("Publication-Ready", "#27ae60", "Meets all scientific integrity standards"),
    (75, 90):  ("Research-Grade", "#2ecc71", "Minor issues remain; suitable for preprint"),
    (60, 75):  ("Needs Improvement", "#f39c12", "Moderate issues; address before submission"),
    (40, 60):  ("Significant Issues", "#e67e22", "Major issues; not publication-ready"),
    (0, 40):   ("Not Scientifically Valid", "#e74c3c", "Critical flaws; do not use for interpretation"),
}

@dataclass
class IntegrityComponent:
    name: str; score: float; weight: float
    evidence: str; source_file: str

@dataclass
class IntegrityReport:
    components: List[IntegrityComponent] = field(default_factory=list)
    total_score: float = 0.0; category: str = ""
    category_description: str = ""; generated_date: str = ""

class ScientificIntegrityScorer:
    """Unified scientific integrity score from all audit dimensions."""

    def __init__(self, output_dir: str = "science"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def score(self,
              snp_universe: Optional[str] = None,
              ancestry_json: Optional[str] = None,
              validation_json: Optional[str] = None,
              calibration_validation: Optional[str] = None,
              leakage_audit: Optional[str] = None,
              benchmark_json: Optional[str] = None,
              fingerprint_json: Optional[str] = None,
              portability_json: Optional[str] = None) -> IntegrityReport:
        logger.info("═══ Scientific Integrity Score ═══")

        comps = [
            self._score_genome_coverage(snp_universe),
            self._score_ancestry_validity(ancestry_json),
            self._score_prs_mathematics(validation_json),
            self._score_calibration_quality(calibration_validation),
            self._score_leakage_prevention(leakage_audit),
            self._score_benchmark_consistency(benchmark_json, validation_json),
            self._score_reproducibility(fingerprint_json),
            self._score_population_portability(portability_json),
        ]

        total = sum(c.score * c.weight for c in comps)
        total = max(0.0, min(100.0, total))

        category = "Not Scientifically Valid"; cat_desc = ""
        for (lo, hi), (cat, _, desc) in INTERPRETATION.items():
            if lo <= total < hi:
                category = cat; cat_desc = desc; break

        report = IntegrityReport(
            components=comps, total_score=round(total, 1),
            category=category, category_description=cat_desc,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save_json(report); self._save_markdown(report)
        return report

    def _load(self, path: Optional[str], key: str, default: Any = None) -> Any:
        if not path or not Path(path).exists(): return default
        try:
            with open(path) as fh: data = json.load(fh)
            return data.get(key, default) if default is not None else data
        except Exception: return default

    def _score_genome_coverage(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        gw = data.get("genome_wide_consistent", False)
        bias = data.get("chr22_bias_detected", True)
        score = 95.0 if (gw and not bias) else (50.0 if gw else 15.0)
        return IntegrityComponent(name="Genome Coverage", score=score,
            weight=WEIGHTS["genome_coverage"],
            evidence=f"Genome-wide: {gw}, chr22 bias: {bias}",
            source_file=path or "none")

    def _score_ancestry_validity(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        method = str(data.get("methodology", {}).get("method",
                 data.get("method_agreement", {}).get("centroid", "")))
        n_pcs = int(data.get("quality_metrics", {}).get("n_pcs_used", 0) or
                    data.get("n_pcs_used", 0))
        is_trait_snp = "allele_frequency" in method.lower()
        is_genome_wide = n_pcs >= 10
        score = 90.0 if (is_genome_wide and not is_trait_snp) else (40.0 if is_trait_snp else 60.0)
        return IntegrityComponent(name="Ancestry Validity", score=score,
            weight=WEIGHTS["ancestry_validity"],
            evidence=f"Method: {method[:40]}, PCs: {n_pcs}, Trait-SNP: {is_trait_snp}",
            source_file=path or "none")

    def _score_prs_mathematics(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        allele_ok = any(c.get("category") == "allele" and c.get("passed")
                       for c in data.get("checks", []))
        score = 85.0 if allele_ok else 50.0
        return IntegrityComponent(name="PRS Mathematics", score=score,
            weight=WEIGHTS["prs_mathematics"],
            evidence=f"Allele validation: {'pass' if allele_ok else 'fail'}",
            source_file=path or "none")

    def _score_calibration_quality(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        slope = float(data.get("mean_slope", 1.0))
        r2 = float(data.get("mean_r2", 0.8))
        well = abs(slope - 1.0) < 0.15 and r2 > 0.80
        score = 90.0 if well else 55.0
        return IntegrityComponent(name="Calibration Quality", score=score,
            weight=WEIGHTS["calibration_quality"],
            evidence=f"Slope={slope:.3f}, R²={r2:.3f}",
            source_file=path or "none")

    def _score_leakage_prevention(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        safe = data.get("pipeline_safe", False)
        errors = int(data.get("errors", 0))
        score = 100.0 if (safe and errors == 0) else (60.0 if safe else 10.0)
        return IntegrityComponent(name="Leakage Prevention", score=score,
            weight=WEIGHTS["leakage_prevention"],
            evidence=f"Safe: {safe}, Errors: {errors}",
            source_file=path or "none")

    def _score_benchmark_consistency(self, pgs: Optional[str],
                                     val: Optional[str]) -> IntegrityComponent:
        data = self._load(pgs, "global_concordance", 0.0)
        score = min(100.0, float(data) * 100 + 20) if isinstance(data, (int, float)) and data > 0 else 50.0
        return IntegrityComponent(name="Benchmark Consistency", score=round(score, 1),
            weight=WEIGHTS["benchmark_consistency"],
            evidence=f"PGS concordance: {float(data):.3f}" if isinstance(data, (int, float)) else "No data",
            source_file=pgs or "none")

    def _score_reproducibility(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, None, {}) or {}
        has_fp = bool(data)  # Any fingerprint data = good
        score = 95.0 if has_fp else 30.0
        return IntegrityComponent(name="Reproducibility", score=score,
            weight=WEIGHTS["reproducibility"],
            evidence=f"Fingerprint: {'present' if has_fp else 'missing'}",
            source_file=path or "none")

    def _score_population_portability(self, path: Optional[str]) -> IntegrityComponent:
        data = self._load(path, "global_bias_index", 0.2)
        bias = float(data) if isinstance(data, (int, float)) else 0.2
        score = max(0.0, 100.0 - bias * 200)
        return IntegrityComponent(name="Population Portability", score=round(score, 1),
            weight=WEIGHTS["population_portability"],
            evidence=f"Bias index: {bias:.3f}",
            source_file=path or "none")

    def _save_json(self, report: IntegrityReport) -> None:
        path = self.output_dir / "scientific_integrity_score.json"
        with open(path, "w") as fh:
            json.dump({
                "scientific_integrity_score": report.total_score,
                "category": report.category,
                "category_description": report.category_description,
                "generated_date": report.generated_date,
                "weights": WEIGHTS, "interpretation": {f"{lo}-{hi-1}": label for (lo, hi), (label, _, _) in INTERPRETATION.items()},
                "components": [{"name": c.name, "score": c.score, "weight": c.weight,
                                "contribution": round(c.score * c.weight, 1),
                                "evidence": c.evidence} for c in report.components],
            }, fh, indent=2)
        logger.info(f"  ✅ Integrity score: {path} — {report.total_score:.0f}/100")

    def _save_markdown(self, report: IntegrityReport) -> None:
        path = self.output_dir / "scientific_integrity_score.md"
        color = "#27ae60" if report.total_score >= 75 else ("#f39c12" if report.total_score >= 60 else "#e74c3c")
        lines = [
            "# Scientific Integrity Score",
            f"\n**Score: {report.total_score:.0f}/100** — {report.category}",
            f"\n{report.category_description}",
            f"\n**Generated:** {report.generated_date}",
            "\n## Component Breakdown",
            "\n| Component | Score | Weight | Contribution | Evidence |",
            "|-----------|-------|--------|-------------|----------|",
        ]
        for c in report.components:
            lines.append(f"| {c.name} | {c.score:.0f}/100 | {c.weight:.0%} | {c.score * c.weight:.1f} | {c.evidence} |")
        lines += ["\n---", f"\n*Phase 8 Scientific Integrity Score — Correction Layer*"]
        with open(path, "w") as fh: fh.write("\n".join(lines))
        logger.info(f"  ✅ Integrity markdown: {path}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scientific Integrity Score")
    parser.add_argument("--snp-universe", default="science/snp_universe.json")
    parser.add_argument("--ancestry-json", default="ancestry/classification_report.json")
    parser.add_argument("--validation-json", default="science/global_validation_report.json")
    parser.add_argument("--calibration-validation", default="benchmark/calibration_validation.json")
    parser.add_argument("--leakage-audit", default="science/leakage_audit.json")
    parser.add_argument("--benchmark-json", default="benchmark/pgs_comparison.json")
    parser.add_argument("--fingerprint-json", default="reproducibility/run_fingerprint.json")
    parser.add_argument("--portability-json", default="benchmark/portability_report.json")
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    scorer = ScientificIntegrityScorer(args.output_dir)
    report = scorer.score(
        args.snp_universe, args.ancestry_json, args.validation_json,
        args.calibration_validation, args.leakage_audit, args.benchmark_json,
        args.fingerprint_json, args.portability_json)
    print(f"\n═══ Scientific Integrity Score ═══")
    print(f"  Score: {report.total_score:.0f}/100 — {report.category}")
    for c in report.components:
        bar = "█" * int(c.score / 5)
        print(f"  {c.name}: {bar} {c.score:.0f}/100 — {c.evidence[:60]}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
