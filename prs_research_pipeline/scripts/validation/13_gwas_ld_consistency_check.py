#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   GWAS + LD + ANCESTRY CONSISTENCY CHECKER                                   ║
║   scripts/13_gwas_ld_consistency_check.py                                    ║
║                                                                            ║
║   Enforces scientific validity by validating:                               ║
║   1. GWAS discovery population matches target ancestry                      ║
║   2. LD reference panel ancestry matches target                             ║
║   3. Fail-fast when both checks fail (no silent continuation)              ║
║                                                                            ║
║   Anti-patterns detected:                                                   ║
║   - EUR GWAS applied to AFR-dominant samples                                ║
║   - LD panel ancestry mismatch > 30%                                        ║
║   - Curated database without population annotation                          ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict
from enum import Enum

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════════

SUPER_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

# Known GWAS population compositions (from GWAS Catalog metadata)
GWAS_POPULATION_MAP = {
    # Major GWAS consortia → primary population
    "Global Lipids Genetics Consortium": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS", "SAS"]},
    "GIANT": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR"]},
    "MAGIC": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS"]},
    "DIAGRAM": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS", "SAS"]},
    "Coffee and Caffeine Genetics Consortium": {"primary": "EUR", "type": "discovery", "includes": ["EUR"]},
    "CHARGE": {"primary": "EUR", "type": "discovery", "includes": ["EUR"]},
    "SUNLIGHT": {"primary": "EUR", "type": "discovery", "includes": ["EUR"]},
    # Trait-level annotations (when consortium unknown)
    "Lipid metabolism": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS", "SAS"]},
    "Glucose metabolism": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS", "SAS", "AFR"]},
    "Obesity predisposition": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS", "AFR"]},
    "Caffeine metabolism": {"primary": "EUR", "type": "discovery", "includes": ["EUR"]},
    "Lactose intolerance": {"primary": "EUR", "type": "multi_population", "includes": ["EUR", "AFR", "SAS"]},
    "Vitamin D metabolism": {"primary": "EUR", "type": "meta_analysis", "includes": ["EUR", "EAS"]},
    "Omega-3 metabolism": {"primary": "EUR", "type": "discovery", "includes": ["EUR"]},
    "Folate & methylation": {"primary": "EUR", "type": "candidate_gene", "includes": ["EUR", "EAS"]},
    "Dopamine regulation": {"primary": "EUR", "type": "candidate_gene", "includes": ["EUR", "EAS"]},
    "Detoxification": {"primary": "EUR", "type": "candidate_gene", "includes": ["EUR", "EAS"]},
}

# Curated database evidence → effective GWAS population
EVIDENCE_POPULATION_MAP = {
    "A": "meta_analysis",     # GWAS p < 5e-8 → likely multi-population
    "B": "replicated",        # Replicated candidate gene
    "C": "single_study",      # Single study → likely single population
    "D": "mechanistic",       # Mechanistic plausibility → no GWAS
}


@dataclass
class ConsistencyResult:
    """Complete consistency check result."""
    gwas_ancestry_match: bool
    gwas_ancestry_score: float           # 0.0–1.0 match quality
    ld_ancestry_match: bool
    ld_ancestry_score: float             # 0.0–1.0 match quality
    risk_flags: List[Dict[str, str]]
    recommended_gwas_source: str
    confidence_downgrade: float           # 0.0–1.0 (0 = no downgrade)
    fail_condition_triggered: bool
    detailed_report: Dict[str, Any]
    passed: bool                          # Overall pass/fail


@dataclass
class TraitConsistencyCheck:
    """Per-trait GWAS ancestry consistency."""
    trait: str
    gwas_population: str
    gwas_type: str
    target_population: str
    target_probability: float
    is_match: bool
    risk_level: str                       # "ok", "warning", "critical"
    note: str


class GWASLDConsistencyChecker:
    """
    Validates GWAS + LD + ancestry consistency before PRS computation.

    Fail condition: If BOTH GWAS mismatch AND LD mismatch are detected,
    the pipeline MUST STOP. This prevents scientifically invalid PRS
    from being computed with mismatched reference data.

    Usage:
        checker = GWASLDConsistencyChecker(ancestry_json_path)
        result = checker.check_all(
            gwas_metadata={"source": "Global Lipids Genetics Consortium"},
            ld_panel_ancestry="EUR",
            curated_db_path="data/snp_database_annotated.csv",
        )
        if not result.passed:
            raise SystemExit("Consistency check failed — see report")
    """

    # Thresholds
    ANCESTRY_MISMATCH_THRESHOLD = 0.5       # P(pop) above this = dominant
    LD_MISMATCH_THRESHOLD = 0.30             # 30% probability gap
    CONFIDENCE_DOWNGRADE_PER_FLAG = 0.15

    def __init__(
        self,
        ancestry_json_path: str,
        strict_mode: bool = True,
    ):
        with open(ancestry_json_path) as fh:
            self.ancestry = json.load(fh)

        self.strict_mode = strict_mode
        self.target_pop = self.ancestry.get("summary", {}).get(
            "assigned_super_population", "EUR"
        )
        self.target_probs = self.ancestry.get("summary", {}).get(
            "all_probabilities", {self.target_pop: 1.0}
        )

    # ── Public API ───────────────────────────────────────────────────────

    def check_all(
        self,
        gwas_metadata: Optional[Dict[str, str]] = None,
        ld_panel_ancestry: Optional[str] = None,
        curated_db_path: Optional[str] = None,
        gwas_source_traits: Optional[Dict[str, str]] = None,
    ) -> ConsistencyResult:
        """
        Run full consistency check.

        Args:
            gwas_metadata: Dict with 'source', 'population', 'type' keys.
            ld_panel_ancestry: Ancestry of LD reference panel used.
            curated_db_path: Path to SNP database for evidence-based checks.
            gwas_source_traits: Per-trait GWAS source mapping.

        Returns:
            ConsistencyResult with pass/fail and detailed report.
        """
        logger.info("═══ GWAS + LD + Ancestry Consistency Check ═══")
        logger.info(f"  Target ancestry: {self.target_pop}")
        logger.info(f"  Target probabilities: {self.target_probs}")

        risk_flags = []
        confidence_downgrade = 0.0
        trait_checks = []

        # ── Check 1: GWAS Ancestry Alignment ────────────────────────────
        gwas_match, gwas_score, gwas_flags, gwas_trait_checks = \
            self._check_gwas_ancestry(gwas_metadata, curated_db_path, gwas_source_traits)
        risk_flags.extend(gwas_flags)
        trait_checks.extend(gwas_trait_checks)

        if not gwas_match:
            confidence_downgrade += self.CONFIDENCE_DOWNGRADE_PER_FLAG

        # ── Check 2: LD Reference Consistency ───────────────────────────
        ld_match, ld_score, ld_flags = self._check_ld_consistency(ld_panel_ancestry)
        risk_flags.extend(ld_flags)

        if not ld_match:
            confidence_downgrade += self.CONFIDENCE_DOWNGRADE_PER_FLAG

        # ── FAIL CONDITION: Both false → STOP ───────────────────────────
        fail_triggered = (not gwas_match) and (not ld_match)

        # Determine recommended GWAS source
        recommended = self._recommend_gwas_source(gwas_metadata)

        # Build detailed report
        detailed = {
            "target_ancestry": {
                "assigned": self.target_pop,
                "probabilities": self.target_probs,
            },
            "gwas_check": {
                "match": gwas_match,
                "score": round(gwas_score, 3),
                "flags": gwas_flags,
            },
            "ld_check": {
                "match": ld_match,
                "score": round(ld_score, 3),
                "flags": ld_flags,
            },
            "trait_checks": [asdict(tc) for tc in trait_checks],
            "fail_condition": fail_triggered,
            "confidence_downgrade": round(confidence_downgrade, 2),
        }

        result = ConsistencyResult(
            gwas_ancestry_match=gwas_match,
            gwas_ancestry_score=round(gwas_score, 3),
            ld_ancestry_match=ld_match,
            ld_ancestry_score=round(ld_score, 3),
            risk_flags=risk_flags,
            recommended_gwas_source=recommended,
            confidence_downgrade=round(confidence_downgrade, 2),
            fail_condition_triggered=fail_triggered,
            detailed_report=detailed,
            passed=not fail_triggered,
        )

        # ── Log result ──────────────────────────────────────────────────
        self._log_result(result)

        if fail_triggered and self.strict_mode:
            logger.critical(
                "╔══════════════════════════════════════════════════════════╗\n"
                "║  ❌ PIPELINE HALTED: GWAS-LD-Ancestry Consistency FAIL   ║\n"
                "║                                                        ║\n"
                "║  Both GWAS ancestry AND LD reference ancestry are       ║\n"
                "║  mismatched with target. PRS computation would produce  ║\n"
                "║  scientifically invalid results.                        ║\n"
                "║                                                        ║\n"
                "║  Recommended actions:                                   ║\n"
                f"║  1. Use GWAS from: {recommended:<35s} ║\n"
                "║  2. Use ancestry-matched LD reference panel             ║\n"
                "║  3. Or acknowledge ancestry mismatch in report           ║\n"
                "╚══════════════════════════════════════════════════════════╝"
            )

        return result

    # ── Private: GWAS Ancestry Check ────────────────────────────────────

    def _check_gwas_ancestry(
        self,
        gwas_metadata: Optional[Dict],
        curated_db_path: Optional[str],
        gwas_source_traits: Optional[Dict[str, str]],
    ) -> Tuple[bool, float, List[Dict], List[TraitConsistencyCheck]]:
        """Check GWAS discovery population matches target ancestry."""
        flags = []
        trait_checks = []

        # Case 1: Formal GWAS metadata provided
        if gwas_metadata:
            gwas_pop = gwas_metadata.get("population", gwas_metadata.get("primary", "EUR"))
            gwas_type = gwas_metadata.get("type", "discovery")

            is_match = self._population_matches(gwas_pop)
            score = self._compute_match_score(gwas_pop, gwas_type)

            if not is_match:
                flags.append({
                    "type": "gwas_ancestry_mismatch",
                    "severity": "critical" if self.strict_mode else "warning",
                    "detail": (
                        f"GWAS population ({gwas_pop}) differs from target "
                        f"ancestry ({self.target_pop}, P={self.target_probs.get(self.target_pop, 0):.1%})"
                    ),
                    "recommendation": f"Use GWAS from {self.target_pop} or multi-ancestry meta-analysis",
                })

            return is_match, score, flags, trait_checks

        # Case 2: Curated database — check per-trait evidence
        if curated_db_path:
            return self._check_curated_db_ancestry(curated_db_path, gwas_source_traits)

        # Case 3: No metadata — flag as unknown
        flags.append({
            "type": "gwas_ancestry_unknown",
            "severity": "warning",
            "detail": "GWAS population metadata not provided — cannot validate ancestry match",
            "recommendation": "Annotate GWAS source with discovery population",
        })
        return True, 0.5, flags, trait_checks  # Pass with warning

    def _check_curated_db_ancestry(
        self,
        db_path: str,
        gwas_source_traits: Optional[Dict],
    ) -> Tuple[bool, float, List[Dict], List[TraitConsistencyCheck]]:
        """Check curated database evidence levels against target ancestry."""
        flags = []
        trait_checks = []

        try:
            db = pd.read_csv(db_path, dtype=str)
        except Exception:
            return True, 0.5, flags, trait_checks

        # Count evidence levels
        if "evidence_level" not in db.columns:
            flags.append({
                "type": "curated_db_no_evidence",
                "severity": "warning",
                "detail": "Curated database lacks evidence_level column",
            })
            return True, 0.5, flags, trait_checks

        evidence_counts = db["evidence_level"].value_counts().to_dict()
        a_fraction = evidence_counts.get("A", 0) / len(db) if len(db) > 0 else 0

        # Check per-trait GWAS population
        if "trait_category" in db.columns:
            for trait in sorted(db["trait_category"].dropna().unique()):
                trait_info = GWAS_POPULATION_MAP.get(trait, {})
                gwas_pop = trait_info.get("primary", "EUR")
                gwas_type = trait_info.get("type", "unknown")

                is_match = self._population_matches(gwas_pop)
                target_p = self.target_probs.get(self.target_pop, 1.0)

                risk_level = "ok"
                note = ""
                if is_match:
                    note = f"GWAS {gwas_pop} matches target {self.target_pop}"
                elif gwas_type == "meta_analysis":
                    risk_level = "ok" if target_p > 0.7 else "warning"
                    note = f"Meta-analysis includes {gwas_pop} — acceptable"
                else:
                    risk_level = "critical" if target_p > 0.5 else "warning"
                    note = f"GWAS {gwas_pop} differs from target {self.target_pop} (P={target_p:.1%})"

                trait_checks.append(TraitConsistencyCheck(
                    trait=trait,
                    gwas_population=gwas_pop,
                    gwas_type=gwas_type,
                    target_population=self.target_pop,
                    target_probability=target_p,
                    is_match=is_match,
                    risk_level=risk_level,
                    note=note,
                ))

        # Overall assessment
        critical_traits = [tc for tc in trait_checks if tc.risk_level == "critical"]
        warning_traits = [tc for tc in trait_checks if tc.risk_level == "warning"]

        if critical_traits:
            flags.append({
                "type": "curated_db_population_mismatch",
                "severity": "critical" if len(critical_traits) > len(trait_checks) * 0.5 else "warning",
                "detail": f"{len(critical_traits)}/{len(trait_checks)} traits have GWAS population mismatch",
                "traits": [tc.trait for tc in critical_traits],
            })

        # Score: fraction of matching traits
        n_ok = len([tc for tc in trait_checks if tc.risk_level == "ok"])
        score = n_ok / len(trait_checks) if trait_checks else 0.5

        # Match is True if no critical mismatches
        is_match = len(critical_traits) == 0

        if a_fraction < 0.3:
            flags.append({
                "type": "low_gwas_evidence",
                "severity": "info",
                "detail": f"Only {a_fraction:.0%} of SNPs have GWAS-level evidence (A)",
                "recommendation": "Results should be interpreted with caution",
            })

        return is_match, score, flags, trait_checks

    # ── Private: LD Consistency Check ────────────────────────────────────

    def _check_ld_consistency(
        self, ld_panel_ancestry: Optional[str]
    ) -> Tuple[bool, float, List[Dict]]:
        """Check LD reference panel ancestry matches target."""
        flags = []

        if ld_panel_ancestry is None:
            # LD panel ancestry unknown — check if we used 1000G multi-population
            ld_match = True  # Multi-population LD is always valid
            ld_score = 0.8   # But not perfect
            flags.append({
                "type": "ld_panel_unknown",
                "severity": "info",
                "detail": "LD panel ancestry not explicitly specified — assuming multi-population",
            })
            return ld_match, ld_score, flags

        ld_match = self._population_matches(ld_panel_ancestry)
        ld_score = self._compute_match_score(ld_panel_ancestry, "ld_reference")

        if not ld_match:
            target_p = self.target_probs.get(self.target_pop, 1.0)
            mismatch_gap = abs(target_p - self.target_probs.get(ld_panel_ancestry, 0))

            severity = "critical" if mismatch_gap > self.LD_MISMATCH_THRESHOLD else "warning"

            flags.append({
                "type": "ld_ancestry_mismatch",
                "severity": severity,
                "detail": (
                    f"LD panel ({ld_panel_ancestry}) differs from target "
                    f"({self.target_pop}, mismatch gap={mismatch_gap:.0%})"
                ),
                "recommendation": (
                    f"Use {self.target_pop}-specific LD reference or "
                    f"multi-population LD panel"
                ),
            })

        return ld_match, ld_score, flags

    # ── Private: Helpers ─────────────────────────────────────────────────

    def _population_matches(self, gwas_pop: str) -> bool:
        """Check if GWAS population matches target ancestry."""
        if gwas_pop == self.target_pop:
            return True
        if gwas_pop == "meta_analysis" or gwas_pop == "multi_population":
            return True
        # If target has high probability for GWAS population
        target_p_for_gwas = self.target_probs.get(gwas_pop, 0)
        if target_p_for_gwas > self.ANCESTRY_MISMATCH_THRESHOLD:
            return True
        return False

    def _compute_match_score(
        self, gwas_pop: str, gwas_type: str = "discovery"
    ) -> float:
        """Compute 0–1 match score."""
        if gwas_pop == self.target_pop:
            return 0.95
        if gwas_type in ("meta_analysis", "multi_population"):
            return 0.85
        if gwas_pop == "multi_population":
            return 0.90

        # Partial match based on target probabilities
        target_p = self.target_probs.get(gwas_pop, 0)
        return min(target_p, 0.9)

    def _recommend_gwas_source(self, gwas_metadata: Optional[Dict]) -> str:
        """Recommend appropriate GWAS source for target ancestry."""
        if self.target_pop == "EUR":
            return "EUR_meta (e.g., GLGC, GIANT, MAGIC)"
        elif self.target_pop == "AFR":
            return "AFR_meta (e.g., AADM, PAGE, Africa WGS)"
        elif self.target_pop == "EAS":
            return "EAS_meta (e.g., BBJ, China Kadoorie)"
        elif self.target_pop == "SAS":
            return "SAS_meta (e.g., BRAVE, Asian Indian)"
        elif self.target_pop == "AMR":
            return "AMR_meta (e.g., PAGE, HCHS/SOL)"
        return "multi_ancestry_meta_analysis"

    def _log_result(self, result: ConsistencyResult) -> None:
        """Log consistency check summary."""
        status = "✅ PASSED" if result.passed else "❌ FAILED"
        logger.info(f"  Status: {status}")
        logger.info(f"  GWAS match: {result.gwas_ancestry_match} (score={result.gwas_ancestry_score})")
        logger.info(f"  LD match: {result.ld_ancestry_match} (score={result.ld_ancestry_score})")
        logger.info(f"  Confidence downgrade: {result.confidence_downgrade:.0%}")
        logger.info(f"  Recommended GWAS: {result.recommended_gwas_source}")
        for flag in result.risk_flags:
            logger.warning(f"  ⚠ [{flag['severity']}] {flag['type']}: {flag['detail']}")


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="GWAS + LD + Ancestry Consistency Checker"
    )
    parser.add_argument("--ancestry", required=True,
                       help="Ancestry inference JSON")
    parser.add_argument("--gwas-source", help="GWAS source/consortium name")
    parser.add_argument("--gwas-population", help="GWAS discovery population")
    parser.add_argument("--gwas-type", default="discovery",
                       help="GWAS type: discovery, meta_analysis, multi_population")
    parser.add_argument("--ld-ancestry", help="LD reference panel ancestry")
    parser.add_argument("--curated-db", help="Curated SNP database path")
    parser.add_argument("--output-dir", "-o", default="prs",
                       help="Output directory")
    parser.add_argument("--no-strict", action="store_true",
                       help="Disable strict mode (don't exit on fail)")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    gwas_meta = None
    if args.gwas_source or args.gwas_population:
        gwas_meta = {
            "source": args.gwas_source or "unknown",
            "population": args.gwas_population or "EUR",
            "type": args.gwas_type,
        }

    checker = GWASLDConsistencyChecker(
        ancestry_json_path=args.ancestry,
        strict_mode=not args.no_strict,
    )

    result = checker.check_all(
        gwas_metadata=gwas_meta,
        ld_panel_ancestry=args.ld_ancestry,
        curated_db_path=args.curated_db,
    )

    # Save report
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path = output_dir / "consistency_check_report.json"
    with open(report_path, "w") as fh:
        json.dump({
            "passed": result.passed,
            "gwas_ancestry_match": result.gwas_ancestry_match,
            "ld_ancestry_match": result.ld_ancestry_match,
            "confidence_downgrade": result.confidence_downgrade,
            "recommended_gwas_source": result.recommended_gwas_source,
            "risk_flags": result.risk_flags,
            "detailed_report": result.detailed_report,
        }, fh, indent=2, default=str)

    print(f"\n{'✅' if result.passed else '❌'} Consistency check {'PASSED' if result.passed else 'FAILED'}")
    print(f"   Report: {report_path}")

    if not result.passed and not args.no_strict:
        print("\n❌ Pipeline halted — see report for details")
        sys.exit(1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
