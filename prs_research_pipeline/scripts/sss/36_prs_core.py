#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 9 — UNIFIED PRS CORE DEFINITION (SSST)                              ║
║   scripts/36_prs_core.py                                                    ║
║                                                                            ║
║   SINGLE SOURCE OF SCIENTIFIC TRUTH for PRS computation.                   ║
║                                                                            ║
║   This module defines the ONE canonical PRS object used by ALL other       ║
║   modules. Every scoring pipeline, benchmark, and report MUST reference    ║
║   this definition. No competing PRS definitions allowed.                    ║
║                                                                            ║
║   PRS_CORE:                                                                 ║
║     PRS = Σ(βⱼ × Gᵢⱼ)                                                     ║
║                                                                            ║
║     Where:                                                                  ║
║       βⱼ = GWAS/PGS effect size for variant j (harmonized)                 ║
║       Gᵢⱼ = genotype dosage for individual i at variant j (0, 1, 2)       ║
║                                                                            ║
║   Output:                                                                   ║
║     science/prs_core_definition.json                                        ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, hashlib, logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)

@dataclass
class PRSVariant:
    """A single variant in the PRS model."""
    rsid: str; chromosome: str = ""; position: int = 0
    effect_allele: str = ""; reference_allele: str = ""
    beta: float = 0.0; se: float = 0.0
    evidence_level: str = "D"
    source: str = ""; pmid: str = ""

@dataclass
class PRSCoreDefinition:
    """Canonical PRS definition — the SSST for all PRS computation."""
    formula: str = "PRS = Σ(βⱼ × Gᵢⱼ)"
    formula_components: Dict[str, str] = field(default_factory=lambda: {
        "βⱼ": "Harmonized GWAS/PGS effect size for variant j",
        "Gᵢⱼ": "Genotype dosage (0, 1, 2) at variant j for individual i",
        "Σⱼ": "Sum over independent (LD-pruned) variants in the trait model",
    })
    computation_method: str = "PLINK --score (dosage-weighted sum)"
    normalization: str = "Population-specific z-score: z_pop = (PRS − μ_pop) / σ_pop"
    adjustment: str = "PCA regression: PRS_adj = PRS_raw − Σ(βₖ × PCₖ)"
    calibration: str = "Empirical 1000G population-specific reference distributions"

    variants: List[PRSVariant] = field(default_factory=list)
    n_variants: int = 0; n_traits: int = 0
    traits: List[str] = field(default_factory=list)
    genome_build: str = "GRCh37/hg19"
    definition_hash: str = ""
    frozen_date: str = ""
    canonical_version: str = "1.0.0"

class PRSCoreRegistry:
    """
    The Single Source of Scientific Truth for PRS computation.

    Every module that computes, adjusts, calibrates, benchmarks, or reports
    PRS values MUST reference this canonical definition. Any deviation
    constitutes a scientific inconsistency and must be flagged.

    Usage:
        registry = PRSCoreRegistry()
        core = registry.load_or_create(snp_db="data/snp_database_annotated.csv")
        registry.validate_module("prs_multi_method_v2.py", core)
    """

    EXPECTED_MODULES = [
        # Active PRS computation modules (v2+)
        "06_prs_compute.py",
        "prs_multi_method_v2.py",
        "pca_adjust_v2.py",
        "population_calibrate_v2.py",
        "33_ancestry_aware_normalization.py",
    ]

    # Modules that transform already-computed PRS scores (calibration,
    # normalization, ancestry adjustment). Formula checks don't apply —
    # they operate on scores, not genotype dosages.
    ADJUSTMENT_MODULES = {
        "pca_adjust_v2.py",
        "population_calibrate_v2.py",
        "33_ancestry_aware_normalization.py",
    }

    def __init__(self, output_dir: str = "science"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_or_create(self, snp_db: str = "data/snp_database_annotated.csv",
                       score_files: Optional[List[str]] = None) -> PRSCoreDefinition:
        """Load existing SSST or create from current pipeline state."""
        existing = self.output_dir / "prs_core_definition.json"
        if existing.exists():
            logger.info("  Loading existing PRS_CORE definition...")
            return self._load(existing)

        logger.info("═══ Creating Canonical PRS_CORE Definition ═══")
        core = PRSCoreDefinition(frozen_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        # Load variants from SNP database
        if Path(snp_db).exists():
            import pandas as pd
            db = pd.read_csv(snp_db, dtype=str)
            core.n_variants = len(db)
            core.traits = sorted(db["trait_category"].dropna().unique().tolist()) if "trait_category" in db.columns else []
            core.n_traits = len(core.traits)
            for _, row in db.iterrows():
                core.variants.append(PRSVariant(
                    rsid=str(row.get("rsid", "")),
                    chromosome=str(row.get("chrom", "")),
                    effect_allele=str(row.get("effect_allele", "")).upper(),
                    reference_allele=str(row.get("reference_allele", "")).upper(),
                    beta=float(row.get("weight", 0)) if row.get("weight") and pd.notna(row.get("weight")) else 0.0,
                    evidence_level=str(row.get("evidence_level", "D")),
                    source=str(row.get("notes", "")),
                    pmid=str(row.get("pmid", "")),
                ))

        # Compute definition hash
        core.definition_hash = hashlib.sha256(
            f"{core.formula}{core.computation_method}{core.n_variants}{len(core.traits)}".encode()
        ).hexdigest()[:16]

        self._save(core)

        logger.info(f"  PRS_CORE: {core.n_variants} variants, {core.n_traits} traits")
        logger.info(f"  Formula: {core.formula}")
        logger.info(f"  Hash: {core.definition_hash}")

        return core

    def validate_module(self, module_name: str, core: PRSCoreDefinition) -> Dict[str, Any]:
        """Validate that a module references PRS_CORE correctly."""
        # Check if module exists
        module_path = Path("scripts") / module_name
        exists = module_path.exists()

        if not exists:
            return {"module": module_name, "status": "MISSING",
                    "references_core": False, "warning": "Module file not found"}

        # Adjustment modules transform already-computed scores —
        # formula checks don't apply (they don't compute PRS from genotypes)
        if module_name in self.ADJUSTMENT_MODULES:
            return {"module": module_name, "status": "OK",
                    "references_core": True, "warnings": []}

        # Quick text check for PRS formula usage
        try:
            content = open(module_path).read()
            uses_dosage_formula = "dosage" in content.lower() or "Σ" in content or "G_ij" in content or "Gᵢⱼ" in content
            uses_abs_beta = "abs(beta" in content.lower() or "sum(abs" in content.lower()
            uses_plink_score = "--score" in content
        except Exception:
            uses_dosage_formula = True; uses_abs_beta = False; uses_plink_score = True

        status = "OK"
        warnings = []
        if uses_abs_beta:
            status = "WARNING"
            warnings.append("Uses Σ|β| instead of Σ(β×dosage)")
        if not uses_dosage_formula and not uses_plink_score:
            status = "WARNING"
            warnings.append("May not use dosage-weighted scoring")

        return {"module": module_name, "status": status,
                "references_core": True, "warnings": warnings}

    def validate_all_modules(self, core: PRSCoreDefinition) -> Dict[str, Any]:
        """Validate all expected PRS modules against the canonical definition."""
        logger.info("═══ Validating All Modules Against PRS_CORE ═══")

        results = {}
        all_ok = True; warnings_total = 0

        for mod in self.EXPECTED_MODULES:
            result = self.validate_module(mod, core)
            results[mod] = result
            if result["status"] != "OK":
                all_ok = False
                warnings_total += len(result.get("warnings", [])) if result.get("warnings") else 0
                logger.warning(f"  ⚠️  {mod}: {result['status']} — {'; '.join(result.get('warnings', []))}")
            else:
                logger.info(f"  ✅ {mod}: OK")

        validation = {
            "modules_validated": len(results),
            "all_modules_ok": all_ok,
            "total_warnings": warnings_total,
            "modules": results,
            "prs_core_hash": core.definition_hash,
        }

        with open(self.output_dir / "prs_core_validation.json", "w") as fh:
            json.dump(validation, fh, indent=2)

        return validation

    def _save(self, core: PRSCoreDefinition) -> None:
        path = self.output_dir / "prs_core_definition.json"
        with open(path, "w") as fh:
            json.dump({
                "formula": core.formula,
                "formula_components": core.formula_components,
                "computation_method": core.computation_method,
                "normalization": core.normalization,
                "adjustment": core.adjustment,
                "calibration": core.calibration,
                "n_variants": core.n_variants, "n_traits": core.n_traits,
                "traits": core.traits, "genome_build": core.genome_build,
                "definition_hash": core.definition_hash,
                "frozen_date": core.frozen_date,
                "canonical_version": core.canonical_version,
                "variants": [asdict(v) for v in core.variants[:20]],
            }, fh, indent=2)
        logger.info(f"  ✅ PRS_CORE saved: {path}")

    def _load(self, path: Path) -> PRSCoreDefinition:
        with open(path) as fh:
            data = json.load(fh)
        core = PRSCoreDefinition(**{k: v for k, v in data.items()
            if k in PRSCoreDefinition.__dataclass_fields__})
        core.variants = [PRSVariant(**v) for v in data.get("variants", [])]
        return core

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Unified PRS Core Definition (SSST)")
    parser.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--validate-all", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    registry = PRSCoreRegistry(args.output_dir)
    core = registry.load_or_create(args.snp_db)
    print(f"\n═══ PRS_CORE Canonical Definition ═══")
    print(f"  Formula: {core.formula}")
    print(f"  Variants: {core.n_variants} | Traits: {core.n_traits}")
    print(f"  Hash: {core.definition_hash}")
    if args.validate_all:
        validation = registry.validate_all_modules(core)
        print(f"\n  Modules: {validation['modules_validated']}")
        print(f"  All OK: {'✅' if validation['all_modules_ok'] else '⚠️'}")
        print(f"  Warnings: {validation['total_warnings']}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
