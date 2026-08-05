#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 9 — UNIFIED PRS RESULT STRUCTURE (SSST)                              ║
║   scripts/37_prs_result_unified.py                                           ║
║                                                                            ║
║   Merges PCA-adjusted, population-calibrated, and ancestry-normalized       ║
║   PRS into a SINGLE canonical output structure.                              ║
║                                                                            ║
║   PRS_RESULT = {                                                             ║
║     raw_score, ancestry_adjusted_score, population_percentile,               ║
║     uncertainty, CI_95, risk_category, calibration_metadata                  ║
║   }                                                                          ║
║                                                                            ║
║   No competing score definitions. All derived views reference this.         ║
║                                                                            ║
║   Output:                                                                    ║
║     prs/PRS_RESULT.json                                                      ║
║     prs/PRS_RESULT.csv                                                       ║
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.constants import PIPELINE_VERSION

logger = logging.getLogger(__name__)

@dataclass
class PRSResultEntry:
    """Single unified PRS result for one trait."""
    trait: str
    raw_score: float
    pca_adjusted_score: float
    ancestry_adjusted_score: float
    population_percentile: float
    population_zscore: float
    uncertainty_score: float
    ci_95_lower: float
    ci_95_upper: float
    risk_category: str  # low, medium, high
    assigned_population: str = "EUR"
    calibration_mu: float = 0.0
    calibration_sigma: float = 1.0
    n_snps_used: int = 0
    n_snps_total: int = 0
    computation_method: str = "PLINK --score"

@dataclass
class UnifiedPRSResult:
    """Complete unified PRS result — the SSST for all downstream consumers."""
    sample_id: str
    pipeline_version: str = PIPELINE_VERSION
    prs_entries: List[PRSResultEntry] = field(default_factory=list)
    ancestry: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    generated_date: str = ""
    result_hash: str = ""

class UnifiedPRSAssembler:
    """
    Assembles all PRS outputs into a single canonical structure.

    Reads from multiple sources (raw, adjusted, calibrated, uncertainty)
    and produces ONE unified PRS_RESULT that all downstream consumers use.

    No module should read individual PRS files after this assembler runs.
    """

    SOURCES = {
        "raw": "prs/prs_raw.csv",
        "pca_adjusted": "prs/pca_adjusted_scores.csv",
        "population_calibrated": "prs/population_calibrated_v2.csv",
        "ancestry_normalized": "prs/ancestry_normalized_scores.csv",
        "uncertainty": "prs/prs_uncertainty.csv",
        "multi_method": "prs/prs_all_methods.csv",
    }

    def __init__(self, output_dir: str = "prs"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def assemble(self, sample_id: str = "SAMPLE_001",
                 ancestry_json: Optional[str] = None,
                 reference_coverage: str = "genome_wide") -> UnifiedPRSResult:
        logger.info("═══ Assembling Unified PRS_RESULT ═══")

        # Load all available sources
        raw = self._load_csv(self.SOURCES["raw"])
        pca_adj = self._load_csv(self.SOURCES["pca_adjusted"])
        pop_cal = self._load_csv(self.SOURCES["population_calibrated"])
        anc_norm = self._load_csv(self.SOURCES["ancestry_normalized"])
        uncert = self._load_csv(self.SOURCES["uncertainty"])
        multi = self._load_csv(self.SOURCES["multi_method"])

        # Load ancestry
        ancestry = {}
        if ancestry_json and Path(ancestry_json).exists():
            with open(ancestry_json) as fh:
                ancestry = json.load(fh)

        assigned_pop = ancestry.get("assigned_population",
            ancestry.get("classification", {}).get("assigned_population", "EUR"))

        # Merge into unified entries
        entries = []
        traits = self._collect_traits(raw, pca_adj, pop_cal)

        for trait in traits:
            entry = PRSResultEntry(trait=trait, raw_score=0, pca_adjusted_score=0, ancestry_adjusted_score=0, population_percentile=50.0, population_zscore=0.0, uncertainty_score=0.0, ci_95_lower=-0.5, ci_95_upper=0.5, risk_category="medium", assigned_population=assigned_pop)

            # Raw score
            raw_row = self._find_trait(raw, trait)
            if raw_row is not None:
                entry.raw_score = round(float(raw_row.get("prs_raw",
                    raw_row.get("prs_adjusted", 0))), 4)
                entry.n_snps_used = int(raw_row.get("n_snps_used", 0))
                entry.n_snps_total = int(raw_row.get("n_snps", 0))

            # PCA-adjusted
            pca_row = self._find_trait(pca_adj, trait)
            if pca_row is not None:
                entry.pca_adjusted_score = round(float(pca_row.get("prs_adjusted",
                    pca_row.get("prs_raw", entry.raw_score))), 4)

            # Population-calibrated
            cal_row = self._find_trait(pop_cal, trait)
            if cal_row is not None:
                entry.population_percentile = round(float(cal_row.get(
                    "percentile_population", 50.0)), 1)
                entry.population_zscore = round(float(cal_row.get(
                    "z_score_population", 0.0)), 4)
                entry.calibration_mu = round(float(cal_row.get("population_mu", 0.0)), 4)
                entry.calibration_sigma = round(float(cal_row.get(
                    "population_sigma", 1.0)), 4)
                entry.risk_category = str(cal_row.get("risk_category", "medium"))

            # Uncertainty
            unc_row = self._find_trait(uncert, trait)
            if unc_row is not None:
                entry.uncertainty_score = round(float(unc_row.get(
                    "uncertainty_score", 0.0)), 4)
                entry.ci_95_lower = round(float(unc_row.get("ci_95_lower", 0.0)), 4)
                entry.ci_95_upper = round(float(unc_row.get("ci_95_upper", 0.0)), 4)

            # Ancestry-adjusted (from ancestry-aware normalization)
            anc_row = self._find_trait(anc_norm, trait)
            if anc_row is not None:
                entry.ancestry_adjusted_score = round(float(anc_row.get(
                    "prs_raw", entry.raw_score)), 4)
                if entry.population_percentile == 0:
                    entry.population_percentile = round(float(anc_row.get(
                        "percentile_population", 50.0)), 1)

            # Default risk if not set
            if not entry.risk_category:
                pctl = entry.population_percentile
                entry.risk_category = "high" if pctl >= 75 else ("low" if pctl <= 25 else "medium")

            entries.append(entry)

        # Build unified result
        result = UnifiedPRSResult(
            sample_id=sample_id,
            prs_entries=entries,
            ancestry={
                "assigned_population": assigned_pop,
                "confidence": ancestry.get("confidence",
                    ancestry.get("classification", {}).get("confidence", "UNKNOWN")),
                "probabilities": ancestry.get("posterior_probabilities",
                    ancestry.get("classification", {}).get("posterior_probabilities", {})),
            },
            metadata={
                "n_traits": len(entries),
                "n_sources_available": sum(1 for s in self.SOURCES.values() if Path(s).exists()),
                "computation_method": "PLINK --score (dosage-weighted)",
                "prs_formula": "PRS = Σ(βⱼ × Gᵢⱼ)",
                "pipeline_version": PIPELINE_VERSION,
                "consolidation_note": "Single Source of Scientific Truth — all PRS values unified",
                "reference_coverage": reference_coverage,
            },
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        )

        result.result_hash = self._hash_result(result)

        self._save_json(result)
        self._save_csv(result)
        self._log_summary(result)

        return result

    def _collect_traits(self, *dataframes) -> List[str]:
        traits = set()
        for df in dataframes:
            if df is not None and len(df) > 0 and "trait" in df.columns:
                traits.update(df["trait"].dropna().unique())
        return sorted(traits)

    def _find_trait(self, df: Optional[pd.DataFrame], trait: str) -> Optional[pd.Series]:
        if df is None or len(df) == 0:
            return None
        if "trait" not in df.columns:
            return df.iloc[0] if len(df) > 0 else None
        matches = df[df["trait"] == trait]
        return matches.iloc[0] if len(matches) > 0 else None

    def _load_csv(self, path: str) -> Optional[pd.DataFrame]:
        if Path(path).exists():
            try:
                return pd.read_csv(path)
            except Exception as e:
                logger.warning(f"  Could not read {path}: {e}")
        return None

    def _hash_result(self, result: UnifiedPRSResult) -> str:
        import hashlib
        sha = hashlib.sha256()
        for e in sorted(result.prs_entries, key=lambda x: x.trait):
            sha.update(f"{e.trait}{e.raw_score}{e.population_percentile}".encode())
        return sha.hexdigest()[:16]

    def _save_json(self, result: UnifiedPRSResult) -> None:
        path = self.output_dir / "PRS_RESULT.json"
        with open(path, "w") as fh:
            json.dump({
                "sample_id": result.sample_id,
                "pipeline_version": result.pipeline_version,
                "generated_date": result.generated_date,
                "result_hash": result.result_hash,
                "ancestry": result.ancestry,
                "metadata": result.metadata,
                "prs_entries": [asdict(e) for e in result.prs_entries],
            }, fh, indent=2)
        logger.info(f"  ✅ PRS_RESULT JSON: {path}")

    def _save_csv(self, result: UnifiedPRSResult) -> None:
        path = self.output_dir / "PRS_RESULT.csv"
        pd.DataFrame([asdict(e) for e in result.prs_entries]).to_csv(path, index=False)
        logger.info(f"  ✅ PRS_RESULT CSV: {path}")

    def _log_summary(self, result: UnifiedPRSResult) -> None:
        logger.info(f"  Entries: {len(result.prs_entries)} traits")
        logger.info(f"  Ancestry: {result.ancestry['assigned_population']}")
        n_high = sum(1 for e in result.prs_entries if e.risk_category == "high")
        n_low = sum(1 for e in result.prs_entries if e.risk_category == "low")
        logger.info(f"  Risk: {n_high} high, {n_low} low")
        for e in result.prs_entries:
            logger.info(f"    {e.trait}: raw={e.raw_score:.3f}, z={e.population_zscore:+.2f}, "
                       f"pctl={e.population_percentile:.0f}%, risk={e.risk_category}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Unified PRS Result Structure")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--ancestry-json", default="ancestry/classification_report.json")
    parser.add_argument("--output-dir", "-o", default="prs")
    parser.add_argument("--reference-coverage", default="genome_wide",
                         choices=["genome_wide", "chr22_only"],
                         help="Whether Stage C ran against the full 1000G reference or fell back to chr22-only")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    assembler = UnifiedPRSAssembler(args.output_dir)
    result = assembler.assemble(args.sample_id, args.ancestry_json, args.reference_coverage)
    print(f"\n═══ Unified PRS_RESULT ═══")
    print(f"  Sample: {result.sample_id}")
    print(f"  Traits: {len(result.prs_entries)}")
    print(f"  Ancestry: {result.ancestry['assigned_population']}")
    print(f"  Hash: {result.result_hash}")
    for e in result.prs_entries[:5]:
        print(f"  {e.trait}: raw={e.raw_score:.3f}, z={e.population_zscore:+.2f}, "
              f"pctl={e.population_percentile:.0f}%, CI=[{e.ci_95_lower:.0f},{e.ci_95_upper:.0f}]")
    return 0

if __name__ == "__main__":
    sys.exit(main())
