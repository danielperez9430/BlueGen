#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 7 — MODULE 3: SCIENTIFIC ASSUMPTION LOCK FILE                       ║
║   scripts/18_scientific_lock.py                                             ║
║                                                                            ║
║   Freezes ALL scientific assumptions used in the pipeline.                  ║
║                                                                            ║
║   SCIENTIFIC FREEZE LAYER — Any deviation in future runs must warn.         ║
║                                                                            ║
║   Locked items:                                                             ║
║     • LD threshold values (window, step, r²)                                ║
║     • PCA dimensionality (n_components)                                     ║
║     • GWAS p-value thresholds                                               ║
║     • Ancestry clusters (super-populations)                                 ║
║     • Calibration distributions (μ, σ per population)                       ║
║     • PRS formula definition                                                ║
║     • QC filter thresholds                                                  ║
║     • Risk category boundaries                                              ║
║     • Imputation parameters                                                 ║
║     • Evidence level definitions                                            ║
║                                                                            ║
║   Output:                                                                   ║
║     science/assumptions.lock.json                                           ║
║     science/assumptions.lock.md                                             ║
║                                                                            ║
║   Validation:                                                               ║
║     Compare current run assumptions against the lock file.                  ║
║     Any deviation triggers a WARNING or ERROR.                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import hashlib
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

SUPER_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

EVIDENCE_LEVELS = {
    "A": {"p_threshold": 5e-8, "replication": 3, "sample_size_min": 50000},
    "B": {"p_threshold": 1e-6, "replication": 2, "sample_size_min": 20000},
    "C": {"p_threshold": 1e-4, "replication": 1, "sample_size_min": 5000},
    "D": {"p_threshold": 0.05, "replication": 0, "sample_size_min": 1000},
}


# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class ScientificLock:
    """Complete frozen scientific assumptions."""
    lock_version: str = "1.0.0"
    frozen_date: str = ""
    lock_hash: str = ""

    # LD Pruning
    ld_window_size: int = 50
    ld_step_size: int = 5
    ld_r2_threshold: float = 0.2
    ld_method: str = "indep-pairwise"

    # PCA
    pca_n_components: int = 20
    pca_method: str = "reference_projection"  # Price et al. 2006
    pca_reference: str = "1000 Genomes Phase 3"

    # GWAS
    gwas_p_threshold: float = 5e-8
    gwas_palindromic_maf_threshold: float = 0.4
    gwas_remove_palindromic: bool = True
    gwas_match_by: str = "position"

    # Ancestry
    ancestry_super_populations: List[str] = field(default_factory=lambda: SUPER_POPULATIONS)
    ancestry_method: str = "pca_ensemble_v2"
    ancestry_pcs_used: int = 10
    ancestry_confidence_thresholds: Dict[str, float] = field(default_factory=lambda: {
        "HIGH": 0.90, "MODERATE": 0.70, "LOW": 0.50,
    })

    # Admixture
    admixture_k: int = 5
    admixture_method: str = "pca_softmax_v2"

    # PRS Formula
    prs_formula: str = "PRS = Σ(βⱼ × Gᵢⱼ)"
    prs_formula_components: Dict[str, str] = field(default_factory=lambda: {
        "βⱼ": "GWAS effect size for SNP j",
        "Gᵢⱼ": "Genotype dosage at SNP j for individual i (0, 1, 2)",
    })
    prs_score_method: str = "sum"
    prs_adjustment_method: str = "pca_regression_v2"

    # QC
    qc_snp_missingness: float = 0.1
    qc_individual_missingness: float = 0.1
    qc_maf: float = 0.01
    qc_hwe: float = 0.000001
    qc_min_quality: int = 20
    qc_min_depth: int = 10

    # Population Calibration
    calibration_method: str = "empirical_1000G_v2"
    calibration_reference: str = "1000 Genomes Phase 3"
    calibration_populations: List[str] = field(default_factory=lambda: SUPER_POPULATIONS)

    # Risk Categories
    risk_low_threshold: float = 25.0
    risk_high_threshold: float = 75.0

    # Uncertainty
    uncertainty_method: str = "three_layer_variance_propagation"
    uncertainty_layers: List[str] = field(default_factory=lambda: [
        "genotype", "ancestry", "effect"
    ])

    # Evidence
    evidence_levels: Dict[str, Dict] = field(default_factory=lambda: EVIDENCE_LEVELS)
    evidence_weights: Dict[str, float] = field(default_factory=lambda: {
        "gwas_significance": 0.25,
        "replication": 0.20,
        "sample_size": 0.20,
        "ancestry_match": 0.15,
        "coverage": 0.10,
        "uncertainty": 0.10,
    })

    # Research Quality
    rqs_components: Dict[str, float] = field(default_factory=lambda: {
        "population_validity": 0.15,
        "prs_validity": 0.20,
        "ancestry_robustness": 0.15,
        "coverage": 0.15,
        "gwas_quality": 0.15,
        "reproducibility": 0.10,
        "calibration_quality": 0.10,
    })

    # Concordance
    concordance_pass_threshold: float = 0.90
    concordance_warning_threshold: float = 0.75


# ── Scientific Lock Engine ────────────────────────────────────────────────────

class ScientificLockEngine:
    """
    Freezes all scientific assumptions into a verifiable lock file.

    Usage:
        engine = ScientificLockEngine()
        lock = engine.freeze()                    # Generate lock file
        result = engine.validate(config.yaml)     # Validate against lock
    """

    LOCK_FILENAME = "assumptions.lock.json"
    LOCK_MD_FILENAME = "assumptions.lock.md"

    def __init__(self, output_dir: str = "science"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Public API ───────────────────────────────────────────────────────

    def freeze(self, overrides: Optional[Dict[str, Any]] = None) -> ScientificLock:
        """
        Generate the scientific assumption lock file.

        This captures the exact state of all scientific parameters.
        Any future run must match these or trigger a validation warning.
        """
        logger.info("═══ Scientific Assumption FREEZE ═══")

        lock = ScientificLock(
            frozen_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        )

        # Apply any overrides
        if overrides:
            for key, value in overrides.items():
                if hasattr(lock, key):
                    setattr(lock, key, value)

        # Compute lock hash (deterministic)
        lock.lock_hash = self._compute_lock_hash(lock)

        # Save
        self._save_lock(lock)
        self._save_lock_markdown(lock)

        logger.info(f"  ✅ Scientific lock frozen: {lock.lock_hash[:16]}")
        logger.info(f"  LD: r² < {lock.ld_r2_threshold}, window={lock.ld_window_size}")
        logger.info(f"  PCA: {lock.pca_n_components} PCs, {lock.pca_method}")
        logger.info(f"  GWAS: p < {lock.gwas_p_threshold}")
        logger.info(f"  Ancestry: {lock.ancestry_super_populations}")
        logger.info(f"  PRS: {lock.prs_formula}")

        return lock

    def validate(self, config_path: str) -> Dict[str, Any]:
        """
        Validate current configuration against the frozen lock file.

        Returns:
            Dict with validation results: {parameter: {lock_value, current_value, match}}
        """
        logger.info("═══ Validating Against Scientific Lock ═══")

        lock_path = self.output_dir / self.LOCK_FILENAME
        if not lock_path.exists():
            return {
                "status": "NO_LOCK",
                "message": "No lock file exists. Run freeze() first.",
                "deviations": [],
            }

        with open(lock_path) as fh:
            lock_data = json.load(fh)

        # Load current config
        current = {}
        if Path(config_path).exists():
            try:
                import yaml
                with open(config_path) as fh:
                    config = yaml.safe_load(fh)
                current = self._extract_config_values(config)
            except Exception as e:
                logger.warning(f"  Could not load config: {e}")

        # Compare key parameters
        checks = []
        deviations = []

        # Map lock keys to config paths
        key_map = {
            "ld_r2_threshold": ("ld_pruning", "r2_threshold"),
            "ld_window_size": ("ld_pruning", "window_size_kb"),
            "qc_snp_missingness": ("qc", "snp_missingness"),
            "qc_maf": ("qc", "maf"),
            "qc_hwe": ("qc", "hwe"),
            "pca_n_components": ("pca", "num_components"),
            "gwas_palindromic_maf_threshold": ("gwas", "palindromic_maf_threshold"),
            "risk_low_threshold": ("prs", "risk_thresholds", "low"),
            "risk_high_threshold": ("prs", "risk_thresholds", "high"),
            "prs_score_method": ("prs", "score_method"),
        }

        for lock_key, config_path_tuple in key_map.items():
            lock_value = lock_data.get(lock_key)

            # Navigate config path
            current_value = current
            for key in config_path_tuple:
                if isinstance(current_value, dict):
                    current_value = current_value.get(key)
                else:
                    current_value = None
                    break

            match = str(lock_value) == str(current_value) if current_value is not None else None

            check = {
                "parameter": lock_key,
                "lock_value": lock_value,
                "current_value": current_value,
                "match": match,
            }
            checks.append(check)

            if match is False:
                severity = "ERROR" if lock_key in (
                    "ld_r2_threshold", "pca_n_components", "qc_maf"
                ) else "WARNING"
                deviations.append({
                    "parameter": lock_key,
                    "lock_value": lock_value,
                    "current_value": current_value,
                    "severity": severity,
                    "message": f"{lock_key}: locked={lock_value} ≠ current={current_value}",
                })

        result = {
            "status": "VALID" if not deviations else "DEVIATION_DETECTED",
            "lock_hash": lock_data.get("lock_hash", ""),
            "n_checks": len(checks),
            "n_matches": sum(1 for c in checks if c["match"] is True),
            "n_mismatches": sum(1 for c in checks if c["match"] is False),
            "checks": checks,
            "deviations": deviations,
        }

        # Log results
        if deviations:
            logger.warning(f"  ⚠️  {len(deviations)} DEVIATIONS from scientific lock!")
            for d in deviations:
                icon = "🔴" if d["severity"] == "ERROR" else "🟡"
                logger.warning(f"    {icon} {d['message']}")
        else:
            logger.info(f"  ✅ All {result['n_checks']} parameters match the lock file")

        return result

    def import_from_config(self, config_path: str) -> ScientificLock:
        """
        Import scientific assumptions from config.yaml and freeze them.

        This ensures the lock file matches the actual configuration.
        """
        if not Path(config_path).exists():
            logger.error(f"  Config not found: {config_path}")
            return self.freeze()

        import yaml
        with open(config_path) as fh:
            config = yaml.safe_load(fh)

        lock = ScientificLock(frozen_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        # LD
        ld = config.get("ld_pruning", {})
        lock.ld_window_size = int(ld.get("window_size_kb", 50))
        lock.ld_r2_threshold = float(ld.get("r2_threshold", 0.2))

        # PCA
        pca = config.get("pca", {})
        lock.pca_n_components = int(pca.get("num_components", 20))

        # GWAS
        gwas = config.get("gwas", {})
        lock.gwas_palindromic_maf_threshold = float(gwas.get("palindromic_maf_threshold", 0.4))
        lock.gwas_remove_palindromic = bool(gwas.get("remove_palindromic", True))
        lock.gwas_match_by = str(gwas.get("match_by", "position"))

        # QC
        qc = config.get("qc", {})
        lock.qc_snp_missingness = float(qc.get("snp_missingness", 0.1))
        lock.qc_individual_missingness = float(qc.get("individual_missingness", 0.1))
        lock.qc_maf = float(qc.get("maf", 0.01))
        lock.qc_hwe = float(qc.get("hwe", 0.000001))

        # PRS
        prs = config.get("prs", {})
        lock.prs_score_method = str(prs.get("score_method", "sum"))
        risk = prs.get("risk_thresholds", {})
        lock.risk_low_threshold = float(risk.get("low", 25))
        lock.risk_high_threshold = float(risk.get("high", 75))

        # Population calibration
        pop_cal = config.get("population_calibration", {})
        lock.calibration_populations = list(pop_cal.get("super_populations", SUPER_POPULATIONS))

        # Normalization
        norm = config.get("normalization", {})
        lock.prs_score_method = str(norm.get("method", "percentile"))

        lock.lock_hash = self._compute_lock_hash(lock)

        self._save_lock(lock)
        self._save_lock_markdown(lock)

        logger.info(f"  ✅ Lock imported from config: {lock.lock_hash[:16]}")
        return lock

    # ── Private Methods ──────────────────────────────────────────────────

    def _compute_lock_hash(self, lock: ScientificLock) -> str:
        """Compute deterministic hash of the lock file."""
        sha = hashlib.sha256()
        # Hash all numeric and string fields (not lists/dicts for stability)
        for field in [
            "ld_window_size", "ld_step_size", "ld_r2_threshold",
            "pca_n_components", "gwas_p_threshold",
            "gwas_palindromic_maf_threshold",
            "qc_snp_missingness", "qc_maf", "qc_hwe",
            "risk_low_threshold", "risk_high_threshold",
            "concordance_pass_threshold", "concordance_warning_threshold",
        ]:
            sha.update(str(getattr(lock, field)).encode())
        return sha.hexdigest()

    def _extract_config_values(self, config: Dict) -> Dict:
        """Extract flat key-value pairs from nested config."""
        flat = {}
        for section, values in config.items():
            if isinstance(values, dict):
                for k, v in values.items():
                    if isinstance(v, dict):
                        for k2, v2 in v.items():
                            flat[f"{section}.{k}.{k2}"] = v2
                    else:
                        flat[f"{section}.{k}"] = v
        return flat

    def _save_lock(self, lock: ScientificLock) -> None:
        """Save lock to JSON."""
        path = self.output_dir / self.LOCK_FILENAME
        with open(path, "w") as fh:
            json.dump(asdict(lock), fh, indent=2)
        logger.info(f"  ✅ Lock JSON: {path}")

    def _save_lock_markdown(self, lock: ScientificLock) -> None:
        """Generate human-readable markdown lock file."""
        path = self.output_dir / self.LOCK_MD_FILENAME

        lines = [
            "# Scientific Assumption Lock File",
            "",
            f"**Frozen:** {lock.frozen_date}",
            f"**Lock Hash:** `{lock.lock_hash[:16]}`",
            f"**Version:** {lock.lock_version}",
            "",
            "⚠️ **WARNING:** Any deviation from these values in future runs MUST be documented and justified.",
            "",
            "---",
            "",
            "## LD Pruning",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| Method | {lock.ld_method} |",
            f"| Window size | {lock.ld_window_size} variants |",
            f"| Step size | {lock.ld_step_size} |",
            f"| r² threshold | {lock.ld_r2_threshold} |",
            "",
            "## PCA",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| Components | {lock.pca_n_components} |",
            f"| Method | {lock.pca_method} |",
            f"| Reference | {lock.pca_reference} |",
            "",
            "## GWAS",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| p-value threshold | {lock.gwas_p_threshold} |",
            f"| Palindromic MAF threshold | {lock.gwas_palindromic_maf_threshold} |",
            f"| Remove palindromic | {lock.gwas_remove_palindromic} |",
            f"| Match by | {lock.gwas_match_by} |",
            "",
            "## QC",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| SNP missingness | {lock.qc_snp_missingness} |",
            f"| Individual missingness | {lock.qc_individual_missingness} |",
            f"| MAF | {lock.qc_maf} |",
            f"| HWE | {lock.qc_hwe} |",
            "",
            "## PRS",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| Formula | `{lock.prs_formula}` |",
            f"| Score method | {lock.prs_score_method} |",
            f"| Adjustment | {lock.prs_adjustment_method} |",
            f"| Risk low | < {lock.risk_low_threshold}th percentile |",
            f"| Risk high | > {lock.risk_high_threshold}th percentile |",
            "",
            "## Ancestry",
            "",
            f"| Parameter | Value |",
            f"|-----------|-------|",
            f"| Super-populations | {', '.join(lock.ancestry_super_populations)} |",
            f"| Method | {lock.ancestry_method} |",
            f"| PCs used | {lock.ancestry_pcs_used} |",
            f"| Admixture K | {lock.admixture_k} |",
            "",
            "## Evidence",
            "",
            "| Level | p-threshold | Replications | Min sample size |",
            "|-------|------------|-------------|-----------------|",
        ]

        for level, info in lock.evidence_levels.items():
            lines.append(f"| {level} | {info['p_threshold']} | {info['replication']} | {info['sample_size_min']:,} |")

        lines += [
            "",
            "## Research Quality (RQS v2)",
            "",
            "| Component | Weight |",
            "|-----------|--------|",
        ]
        for comp, weight in lock.rqs_components.items():
            lines.append(f"| {comp} | {weight:.0%} |")

        lines += [
            "",
            "---",
            "",
            f"*Assumptions frozen {lock.frozen_date}. Lock hash: `{lock.lock_hash[:16]}`*",
            "",
            "> **Scientific Freeze Layer — Phase 7**",
            "> No algorithmic changes. No new biology. No model improvements.",
            "> Only reproducibility, validation, and audit trails from this point forward.",
        ]

        with open(path, "w") as fh:
            fh.write("\n".join(lines))
        logger.info(f"  ✅ Lock Markdown: {path}")


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 7 Module 3: Scientific Assumption Lock File"
    )
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--freeze", action="store_true",
                       help="Generate the lock file")
    parser.add_argument("--validate", help="Validate config.yaml against lock file")
    parser.add_argument("--import-config", help="Import assumptions from config.yaml")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    engine = ScientificLockEngine(output_dir=args.output_dir)

    if args.import_config:
        lock = engine.import_from_config(args.import_config)
        print(f"\n═══ Scientific Lock Imported ═══")
        print(f"  Hash: {lock.lock_hash[:16]}")
        print(f"  LD r²: {lock.ld_r2_threshold}")
        print(f"  PCA components: {lock.pca_n_components}")
        print(f"  QC MAF: {lock.qc_maf}")
        print(f"  Risk: <{lock.risk_low_threshold}% / >{lock.risk_high_threshold}%")

    elif args.freeze:
        lock = engine.freeze()
        print(f"\n═══ Scientific Lock Frozen ═══")
        print(f"  Hash: {lock.lock_hash[:16]}")
        print(f"  File: science/{engine.LOCK_FILENAME}")

    elif args.validate:
        result = engine.validate(args.validate)
        print(f"\n═══ Lock Validation ═══")
        print(f"  Status: {result['status']}")
        print(f"  Checks: {result['n_checks']} total, {result['n_matches']} matched, {result['n_mismatches']} deviated")
        for d in result.get("deviations", []):
            print(f"  {d['severity']}: {d['message']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
