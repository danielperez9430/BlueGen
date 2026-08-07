"""
Pipeline data loading extracted from comprehensive_report.py::main()
(IMPROVEMENT_PLAN.md 1.6, Phase 1) - moved verbatim, no logic changes.
"""

import json
import logging
import os
import sys
from pathlib import Path

import pandas as pd

from report.interpretations import load_json, safe_float

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bluegen.schemas import (  # noqa: E402
    PRSResult, AncestryModel, CalibrationValidationReport, UncertaintyReport,
    load_and_validate_json,
)

logger = logging.getLogger(__name__)


def _load_validated(path: str, model_cls) -> dict:
    """Read-time validation for one of the 4 pydantic-schema'd critical
    artifacts (IMPROVEMENT_PLAN.md 2.2). A missing file is a legitimate
    partial-pipeline state (same as load_json's default={}), not a shape
    bug, so it's not an error here. A file that exists but doesn't match
    its producer's real contract raises pydantic.ValidationError, failing
    loudly at this single choke point instead of every downstream
    .get(key, default) silently degrading.
    """
    if not os.path.exists(path):
        return {}
    return load_and_validate_json(model_cls, path).model_dump(mode="json")


def load_report_data(args) -> dict:
    """Load every pipeline JSON/CSV output `build_html_report()` needs.

    `args` is the argparse.Namespace from comprehensive_report.py::main()
    (needs .snp_db at minimum). Returns the same `data` dict shape main()
    has always built, including the synthetic `_cal_lookup`/`_uncert_lookup`/
    `_evidence_lookup`/`_recommendation_lookup`/`_polarity_inverted`/
    `_pgs_coverage_lookup` keys.
    """
    # Load all data sources
    logger.info("═══ Loading Pipeline Data ═══")

    data = {}

    # Core outputs
    data["prs_result"] = _load_validated("prs/PRS_RESULT.json", PRSResult)
    data["ancestry"] = _load_validated("science/ANCESTRY_MODEL.json", AncestryModel)
    data["validation"] = load_json("science/global_validation_report.json")
    data["integrity"] = load_json("FINAL_SCIENTIFIC_SCORE.json")
    data["benchmark"] = load_json("benchmark/VALIDATION_REPORT.json")
    data["adversarial"] = load_json("science/adversarial_validation_report.json")
    data["failure_map"] = load_json("science/failure_mode_map.json")
    data["leakage"] = load_json("science/pipeline_gate_check.json")
    data["quality_delta"] = load_json("benchmark/quality_delta.json")
    # Extended data sources
    data["uncertainty_report"] = _load_validated("prs/uncertainty_report.json", UncertaintyReport)
    # prs/population_calibration_report.json has no producer anywhere in the
    # codebase (confirmed via git log -S across all history) - it's a
    # phantom artifact, not schema'd, read tolerantly like everything else
    # below (IMPROVEMENT_PLAN.md 2.2).
    data["calibration_report"] = load_json("prs/population_calibration_report.json")
    data["gwas_consortium"] = load_json("benchmark/gwas_consortium_validation.json")
    data["portability"] = load_json("benchmark/portability_report.json")
    data["reproducibility"] = load_json("reproducibility/run_fingerprint.json")
    data["consistency"] = load_json("prs/consistency_check_report.json")
    data["leakage_audit"] = load_json("science/leakage_audit.json")
    data["snp_universe"] = load_json("science/snp_universe.json")
    data["pgs_calibration"] = load_json("prs/pgs_scores/pgs_calibration_report.json")
    data["clinvar"] = load_json("clinvar/clinvar_pathogenic_variants.json")
    data["pharmgkb"] = load_json("pharmgkb/pharmgkb_drug_report.json")
    data["deep_ancestry"] = load_json("ancestry/deep_ancestry.json")

    # NEW: Load subcontinental assignment if available
    subc_data = load_json("pca/subcontinental_assignment.json")
    if subc_data and subc_data.get("assigned_sub_population"):
        # Merge into deep_ancestry for display
        if data["deep_ancestry"]:
            data["deep_ancestry"]["sub_continental"] = subc_data
        else:
            data["deep_ancestry"] = {"sub_continental": subc_data}
        logger.info(f"  ✅ Subcontinental assignment: {subc_data.get('assigned_sub_population')}")

    # NEW: Load calibration validation data per trait (for confidence scores)
    data["calibration_validation"] = _load_validated(
        "benchmark/calibration_validation.json", CalibrationValidationReport)

    # Build per-trait calibration lookup
    cal_lookup = {}
    for v in data["calibration_validation"].get("validations", []):
        cal_lookup[v["trait"].lower()] = v
    data["_cal_lookup"] = cal_lookup

    # Build per-trait uncertainty decomposition lookup
    uncert_lookup = {}
    for r in data.get("uncertainty_report", {}).get("results", []):
        uncert_lookup[r["trait"].lower()] = r
    data["_uncert_lookup"] = uncert_lookup

    # Build per-trait evidence level lookup from SNP database
    evidence_lookup = {}
    snp_db_path = args.snp_db
    if os.path.exists(snp_db_path):
        try:
            snp_db = pd.read_csv(snp_db_path, dtype=str)
            for col in ["trait_category", "trait", "Trait", "trait_name"]:
                if col in snp_db.columns:
                    trait_col = col
                    break
            else:
                trait_col = None
            if trait_col:
                ev_map = {"A": 100, "B": 75, "C": 50, "D": 25}
                for _, row in snp_db.iterrows():
                    trait = str(row.get(trait_col, "")).strip().lower()
                    ev = str(row.get("evidence_level", row.get("evidence", "C"))).strip().upper()
                    ev_score = ev_map.get(ev, 50)
                    if trait not in evidence_lookup:
                        evidence_lookup[trait] = []
                    evidence_lookup[trait].append(ev_score)
        except Exception:
            pass
    data["_evidence_lookup"] = evidence_lookup

    # Curated, evidence-cited per-trait action recommendations (IMPROVEMENT_PLAN.md
    # 1.4). Deliberately a small hand-verified subset, not auto-generated - see
    # data/trait_recommendations.json's _meta for the curation scope/rationale.
    #
    # _polarity_inverted lists traits where risk_category=="high" means a
    # favorable result, not elevated risk (e.g. Morning chronotype, Cognitive
    # function) - see risk_badge()'s docstring. Both this and
    # recommendation_lookup are keyed by lowercased trait name; any other
    # top-level "_"-prefixed key (just _meta today) is metadata, not a trait.
    recommendation_lookup = {}
    polarity_inverted = set()
    trait_reco_path = os.path.join(os.path.dirname(args.snp_db), "trait_recommendations.json")
    if os.path.exists(trait_reco_path):
        try:
            with open(trait_reco_path) as fh:
                raw = json.load(fh)
            recommendation_lookup = {k.lower(): v for k, v in raw.items() if not k.startswith("_")}
            polarity_inverted = {t.lower() for t in raw.get("_polarity_inverted", [])}
        except Exception:
            pass
    data["_recommendation_lookup"] = recommendation_lookup
    data["_polarity_inverted"] = polarity_inverted

    # Build PGS coverage lookup from pgs_results.csv
    pgs_coverage_lookup = {}
    pgs_results_path = "prs/pgs_scores/pgs_results.csv"
    if os.path.exists(pgs_results_path):
        try:
            pgs_results = pd.read_csv(pgs_results_path, dtype=str)
            for _, row in pgs_results.iterrows():
                pgs_id = str(row.get("pgs_id", "")).strip()
                n_used = safe_float(row.get("n_snps_used", 0))
                n_total = safe_float(row.get("n_snps_in_score", 0))
                if pgs_id:
                    pgs_coverage_lookup[pgs_id] = {"n_used": n_used, "n_total": n_total}
        except Exception:
            pass
    data["_pgs_coverage_lookup"] = pgs_coverage_lookup

    # Log what was found
    for name, d in data.items():
        status = "✅" if d else "⬚ (missing)"
        logger.info(f"  {status} {name}")

    # Fall back to ancestry inference if ANCESTRY_MODEL is empty or UNKNOWN
    anc_pop = data["ancestry"].get("assigned_population", "UNKNOWN")
    if not data["ancestry"] or anc_pop in ("UNKNOWN", None, ""):
        alt = load_json("pca/ancestry_inference.json")
        if alt:
            # Map ancestry_inference format to ANCESTRY_MODEL format
            summary = alt.get("summary", {})
            data["ancestry"] = {
                "assigned_population": summary.get("assigned_super_population", "EUR"),
                "posterior_probabilities": summary.get("all_probabilities", {}),
                "confidence": summary.get("confidence", "MODERATE"),
                "n_reference_samples": alt.get("methodology", {}).get("snps_used", 2504),
                "n_pcs": 20,
                "method": alt.get("methodology", {}).get("method", "allele_frequency_distance"),
            }
            logger.info(f"  Using pca/ancestry_inference.json: {data['ancestry']['assigned_population']} ({data['ancestry']['confidence']})")

    return data
