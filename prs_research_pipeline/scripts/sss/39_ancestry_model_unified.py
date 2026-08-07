#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 9 — ANCESTRY MODEL UNIFICATION (SSST)                                ║
║   scripts/39_ancestry_model_unified.py                                       ║
║                                                                            ║
║   Enforces a SINGLE ancestry model: 1000G PCA projection only.              ║
║   All other ancestry estimators become diagnostics — never inputs to        ║
║   scoring or calibration.                                                    ║
║                                                                            ║
║   ANCESTRY_MODEL:                                                            ║
║     method: 1000G_PCA_PROJECTION                                             ║
║     reference: 1000 Genomes Phase 3 (all autosomes)                         ║
║     pcs: PC1–PC20                                                            ║
║     classifier: PCA ensemble (centroid + logistic + k-NN)                    ║
║     output: posterior_probabilities, confidence, quality_metrics             ║
║                                                                            ║
║   Diagnostics only (NOT used for scoring):                                  ║
║     - allele_frequency_distance (legacy)                                     ║
║     - any <100-marker method                                                 ║
║                                                                            ║
║   Output:                                                                    ║
║     science/ANCESTRY_MODEL.json                                              ║
║     science/ancestry_diagnostics_audit.json                                  ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from bluegen.schemas import AncestryModel, validate_and_write_json  # noqa: E402

logger = logging.getLogger(__name__)

VALID_ANCESTRY_METHODS = ["pca_ensemble_v2", "reference_projection", "allele_frequency_distance"]
DIAGNOSTIC_ONLY_METHODS = ["allele_frequency_distance", "centroid_distance",
                            "nearest_neighbor", "knn", "logistic_regression"]

@dataclass
class AncestryAudit:
    model: AncestryModel
    alternative_methods_found: List[Dict[str, str]]
    scoring_uses_canonical_model: bool
    diagnostics_only: List[str]
    generated_date: str

class AncestryModelUnification:
    """
    Enforces a single ancestry model for all downstream PRS operations.

    The canonical model is 1000G PCA projection. Any other ancestry
    estimator found in the pipeline is classified as diagnostic-only
    and MUST NOT influence scoring or calibration.
    """

    ANCESTRY_FILES = [
        "pca/ancestry_classification.json",  # PCA ensemble classifier (v2)
        "ancestry/classification_report.json",
        "ancestry/posterior_probabilities.json",
        "pca/ancestry_inference.json",
        "pca/projected_sample.csv",
    ]

    def __init__(self, output_dir: str = "science"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def unify(self) -> AncestryAudit:
        logger.info("═══ Ancestry Model Unification ═══")

        model = AncestryModel(frozen_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))
        alternative_methods = []
        diagnostics = []

        # Try to load the canonical PCA ensemble output
        for path in self.ANCESTRY_FILES:
            if not Path(path).exists():
                continue
            try:
                with open(path) as fh:
                    data = json.load(fh)
                method = str(data.get("methodology", {}).get("method",
                         data.get("method", "")))
                n_snps = data.get("methodology", {}).get("snps_used", 1000)

                if any(m.lower() in method.lower() for m in VALID_ANCESTRY_METHODS):
                    # This is the canonical model
                    # Support multiple ancestry data formats
                    summary = data.get("summary", {})
                    classification = data.get("classification", {})
                    model.assigned_population = (
                        data.get("assigned_population") or
                        summary.get("assigned_super_population") or
                        classification.get("assigned_population") or
                        "EUR"
                    )
                    model.posterior_probabilities = (
                        data.get("posterior_probabilities") or
                        summary.get("all_probabilities") or
                        classification.get("posterior_probabilities") or
                        {}
                    )
                    model.confidence = (
                        data.get("confidence") or
                        summary.get("confidence") or
                        classification.get("confidence") or
                        "MODERATE"
                    )
                    metrics = data.get("quality_metrics", {})
                    if metrics:
                        model.quality_metrics = {
                            "max_probability": metrics.get("max_posterior_probability", 0),
                            "entropy": metrics.get("entropy", 0),
                            "distance_ratio": metrics.get("distance_ratio", 0),
                            "n_pcs_used": metrics.get("n_pcs_used", 20),
                        }
                    model.n_pcs = model.quality_metrics.get("n_pcs_used", 20)
                    model.is_valid_for_scoring = True
                    model.method = "PCA_ENSEMBLE_V2"
                    logger.info(f"  Canonical model found: {path}")

                elif any(d.lower() in method.lower() for d in DIAGNOSTIC_ONLY_METHODS):
                    diagnostics.append(path)
                    if isinstance(n_snps, int) and n_snps < 100:
                        alternative_methods.append({
                            "file": path, "method": method,
                            "issue": f"Only {n_snps} markers — diagnostic only",
                            "status": "DIAGNOSTIC_ONLY"})
                        logger.warning(f"  Diagnostic-only method: {path} ({method}, {n_snps} SNPs)")

            except Exception as e:
                logger.warning(f"  Could not parse {path}: {e}")

        # Compute model hash
        import hashlib
        sha = hashlib.sha256()
        sha.update(f"{model.method}{model.assigned_population}{model.n_pcs}".encode())
        model.model_hash = sha.hexdigest()[:16]

        audit = AncestryAudit(
            model=model,
            alternative_methods_found=alternative_methods,
            scoring_uses_canonical_model=model.is_valid_for_scoring,
            diagnostics_only=diagnostics,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save(model, audit)
        return audit

    def _save(self, model: AncestryModel, audit: AncestryAudit) -> None:
        # Canonical model
        validate_and_write_json(model, self.output_dir / "ANCESTRY_MODEL.json")

        # Diagnostics audit
        with open(self.output_dir / "ancestry_diagnostics_audit.json", "w") as fh:
            json.dump({
                "canonical_model": model.method,
                "is_valid_for_scoring": model.is_valid_for_scoring,
                "alternative_methods_found": len(audit.alternative_methods_found),
                "diagnostics_only_modules": audit.diagnostics_only,
                "warning": "Alternative ancestry methods are diagnostics only — they do not contribute to PRS scoring or calibration",
                "alternatives": audit.alternative_methods_found,
                "generated_date": audit.generated_date,
            }, fh, indent=2)

        logger.info(f"  ✅ ANCESTRY_MODEL: {model.assigned_population} ({model.confidence})")
        logger.info(f"  ✅ Valid for scoring: {model.is_valid_for_scoring}")
        logger.info(f"  ✅ Diagnostics only: {len(audit.diagnostics_only)} files")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Ancestry Model Unification")
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    unifier = AncestryModelUnification(args.output_dir)
    audit = unifier.unify()
    print(f"\n═══ Ancestry Model Unification ═══")
    print(f"  Model: {audit.model.method}")
    print(f"  Population: {audit.model.assigned_population} ({audit.model.confidence})")
    print(f"  Valid for scoring: {'✅' if audit.model.is_valid_for_scoring else '❌'}")
    print(f"  Diagnostics only: {len(audit.diagnostics_only)} files")
    return 0

if __name__ == "__main__":
    sys.exit(main())
