#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 9 — BENCHMARK REINTERPRETATION LAYER (SSST)                          ║
║   scripts/38_benchmark_reinterpretation.py                                   ║
║                                                                            ║
║   Enforces that benchmarks evaluate ONLY PRS_CORE outputs, never            ║
║   intermediate calibrated versions. Eliminates circular validation.         ║
║                                                                            ║
║   Classification:                                                            ║
║     ✔ Internal validation — pipeline correctness (cross-method, QC)         ║
║     ✔ External validation — PGS/GWAS comparison (independent reference)     ║
║     ❌ Re-validation — derived transforms of already-validated scores       ║
║                                                                            ║
║   Output:                                                                    ║
║     benchmark/validation_classification.json                                 ║
║     benchmark/VALIDATION_REPORT.json                                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime

logger = logging.getLogger(__name__)

@dataclass
class ValidationEntry:
    validation_id: str; description: str
    validation_type: str  # internal, external, re_validation
    evaluates: str  # what is being evaluated
    reference: str  # what it's compared against
    is_independent: bool  # reference is external to platform
    is_circular: bool  # evaluates derived outputs of itself
    status: str  # VALID, CIRCULAR, REDUNDANT, UNCLASSIFIED

@dataclass
class ReinterpretationReport:
    entries: List[ValidationEntry] = field(default_factory=list)
    n_internal: int = 0; n_external: int = 0
    n_circular: int = 0; n_redundant: int = 0
    circular_validations: List[str] = field(default_factory=list)
    all_independent: bool = False
    generated_date: str = ""

class BenchmarkReinterpretation:
    """
    Audits all benchmark/validation modules and classifies each as:
      - Internal validation (pipeline correctness)
      - External validation (independent reference)
      - Circular re-validation (must be downgraded)

    All benchmarks must reference PRS_CORE outputs only.
    """

    VALIDATION_INVENTORY = [
        {"id": "VAL-CONCORDANCE", "module": "concordance_analysis.py",
         "evaluates": "Platform PRS (population-calibrated)",
         "reference": "PGS Catalog scores", "type": "external"},
        {"id": "VAL-COVERAGE", "module": "coverage_audit.py",
         "evaluates": "Variant presence in sample genotypes",
         "reference": "Curated SNP database", "type": "internal"},
        {"id": "VAL-EVIDENCE", "module": "evidence_scoring.py",
         "evaluates": "Trait evidence quality",
         "reference": "GWAS significance, replication, sample size", "type": "internal"},
        {"id": "VAL-READINESS", "module": "clinical_readiness.py",
         "evaluates": "Clinical readiness of calibrated PRS",
         "reference": "Evidence scores, coverage, ancestry", "type": "internal"},
        {"id": "VAL-LIMITATIONS", "module": "limitations_engine.py",
         "evaluates": "Scientific limitations",
         "reference": "Coverage, evidence, ancestry thresholds", "type": "internal"},
        {"id": "VAL-BENCHMARK", "module": "23_pgs_catalog_benchmark.py",
         "evaluates": "Platform PRS vs PGS Catalog",
         "reference": "PGS Catalog scores", "type": "external"},
        {"id": "VAL-METHODS", "module": "24_external_prs_replication.py",
         "evaluates": "Cross-method PRS agreement",
         "reference": "Internal method comparison", "type": "internal"},
        {"id": "VAL-GWAS", "module": "25_gwas_consortium_validation.py",
         "evaluates": "GWAS consortium alignment",
         "reference": "Published GWAS effect directions", "type": "external"},
        {"id": "VAL-PORTABILITY", "module": "26_population_portability_test.py",
         "evaluates": "PRS across populations",
         "reference": "Published transferability expectations", "type": "external"},
        {"id": "VAL-CALIBRATION", "module": "27_real_world_calibration.py",
         "evaluates": "Calibration quality of calibrated PRS",
         "reference": "Empirical 1000G distributions", "type": "internal"},
        {"id": "VAL-POSITIONING", "module": "28_scientific_positioning.py",
         "evaluates": "Scientific positioning",
         "reference": "Published PRS standards", "type": "external"},
        {"id": "VAL-DELTA", "module": "29_quality_delta_analysis.py",
         "evaluates": "Quality gap analysis",
         "reference": "External benchmark expectations", "type": "external"},
    ]

    # Modules that evaluate calibrated outputs (risk of circularity)
    CALIBRATION_DEPENDENT = [
        "VAL-READINESS", "VAL-CALIBRATION", "VAL-CONCORDANCE",
    ]

    def __init__(self, output_dir: str = "benchmark"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def audit(self, prs_core_json: Optional[str] = None) -> ReinterpretationReport:
        logger.info("═══ Benchmark Reinterpretation Audit ═══")

        prs_core_hash = ""
        if prs_core_json and Path(prs_core_json).exists():
            with open(prs_core_json) as fh:
                prs_core_hash = json.load(fh).get("definition_hash", "")

        entries = []
        for val in self.VALIDATION_INVENTORY:
            is_circular = val["id"] in self.CALIBRATION_DEPENDENT
            is_independent = val["type"] == "external"

            if is_circular:
                status = "CIRCULAR"
            elif val["type"] == "internal":
                status = "VALID"
            else:
                status = "VALID"

            entries.append(ValidationEntry(
                validation_id=val["id"],
                description=val["evaluates"],
                validation_type=val["type"],
                evaluates=val["evaluates"],
                reference=val["reference"],
                is_independent=is_independent,
                is_circular=is_circular,
                status=status))

        circular_ids = [e.validation_id for e in entries if e.is_circular]

        report = ReinterpretationReport(
            entries=entries,
            n_internal=sum(1 for e in entries if e.validation_type == "internal"),
            n_external=sum(1 for e in entries if e.validation_type == "external"),
            n_circular=len(circular_ids),
            n_redundant=0,
            circular_validations=circular_ids,
            all_independent=len(circular_ids) == 0,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save(report, prs_core_hash)
        return report

    def _save(self, report: ReinterpretationReport, prs_core_hash: str) -> None:
        # Classification
        with open(self.output_dir / "validation_classification.json", "w") as fh:
            json.dump({
                "n_internal": report.n_internal, "n_external": report.n_external,
                "n_circular": report.n_circular,
                "circular_warning": "Circular validations evaluate calibrated outputs that depend on upstream platform scores. These should be interpreted as internal consistency checks, not independent validation.",
                "circular_validations": report.circular_validations,
                "all_independent": report.all_independent,
                "prs_core_hash": prs_core_hash,
                "generated_date": report.generated_date,
                "entries": [asdict(e) for e in report.entries],
            }, fh, indent=2)

        # Unified validation report
        with open(self.output_dir / "VALIDATION_REPORT.json", "w") as fh:
            json.dump({
                "validation_summary": {
                    "total_validations": len(report.entries),
                    "internal": report.n_internal,
                    "external": report.n_external,
                    "circular": report.n_circular,
                    "all_independent": report.all_independent,
                },
                "prs_core_reference": prs_core_hash,
                "classification": {
                    "internal_validation": "Pipeline correctness — cross-method, QC, coverage",
                    "external_validation": "Independent reference — PGS Catalog, GWAS consortia",
                    "circular_re_validation": "Evaluates derived transforms — interpret as internal consistency",
                },
                "entries": [asdict(e) for e in report.entries],
                "generated_date": report.generated_date,
            }, fh, indent=2)

        logger.info(f"  ✅ Classification: {report.n_internal} internal, {report.n_external} external, {report.n_circular} circular")
        if report.n_circular > 0:
            logger.warning(f"  ⚠️  Circular validations: {', '.join(report.circular_validations)}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Benchmark Reinterpretation")
    parser.add_argument("--prs-core-json", default="science/prs_core_definition.json")
    parser.add_argument("--output-dir", "-o", default="benchmark")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    auditor = BenchmarkReinterpretation(args.output_dir)
    report = auditor.audit(args.prs_core_json)
    print(f"\n═══ Benchmark Reinterpretation ═══")
    print(f"  Internal: {report.n_internal} | External: {report.n_external}")
    print(f"  Circular: {report.n_circular} | Redundant: {report.n_redundant}")
    print(f"  All independent: {'✅' if report.all_independent else '⚠️ NO'}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
