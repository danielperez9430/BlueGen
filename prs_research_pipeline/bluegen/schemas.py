"""
Pydantic schemas for critical pipeline JSON artifacts (IMPROVEMENT_PLAN.md 2.2).

These are shape contracts, not business-rule validators: no range/bounds
checks (e.g. percentile clamping) - only "does this have the fields real
producers already emit, with the right types." Fields are required only
where every real producer always sets them and no consumer defensively
omits them (confirmed via a full-repo grep of every reader of these four
files - the only strict `dict["key"]` site across all of them is `trait`
on each array entry); everything else is optional with the same default
the dataclass-based producers already used, so validating today's real
output never rejects it.

Unknown fields are preserved (extra="allow") rather than rejected or
silently dropped - a producer adding a field later shouldn't need a
schema change to keep working, and a consumer relying on model_dump()
shouldn't lose data it didn't know to ask for.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple, Type, TypeVar, Union

from pydantic import BaseModel, ConfigDict

ModelT = TypeVar("ModelT", bound=BaseModel)


class _LooseModel(BaseModel):
    model_config = ConfigDict(extra="allow")


# ─────────────────────────── PRS_RESULT.json ────────────────────────────
# Producer: scripts/sss/37_prs_result_unified.py

class PRSResultEntry(_LooseModel):
    trait: str
    raw_score: float = 0.0
    pca_adjusted_score: float = 0.0
    ancestry_adjusted_score: float = 0.0
    population_percentile: float = 50.0
    population_zscore: float = 0.0
    uncertainty_score: float = 0.0
    ci_95_lower: float = -0.5
    ci_95_upper: float = 0.5
    risk_category: str = "medium"
    assigned_population: str = "EUR"
    calibration_mu: float = 0.0
    calibration_sigma: float = 1.0
    n_snps_used: int = 0
    n_snps_total: int = 0
    computation_method: str = "PLINK --score"


class PRSResult(_LooseModel):
    sample_id: str
    pipeline_version: str = ""
    generated_date: str = ""
    result_hash: str = ""
    ancestry: Dict[str, Any] = {}
    metadata: Dict[str, Any] = {}
    prs_entries: List[PRSResultEntry] = []


# ─────────────────────────── ANCESTRY_MODEL.json ─────────────────────────
# Producers: scripts/sss/39_ancestry_model_unified.py (Producer A) AND
#            scripts/stages/pca_ancestry_classifier.py (Producer B).
# quality_metrics stays Dict[str, Any] on purpose: the two producers write
# structurally different sub-shapes (A: max_probability/entropy/
# distance_ratio/n_pcs_used; B: confidence_ratio/knn_mean_distance) and
# forcing one shape on both would invent a contract neither actually
# follows.

class AncestryModel(_LooseModel):
    method: str = "PCA_ENSEMBLE_V2"
    reference_panel: str = "1000 Genomes Phase 3 (all autosomes)"
    n_pcs: int = 20
    n_reference_samples: int = 2504
    super_populations: List[str] = ["EUR", "AFR", "EAS", "SAS", "AMR"]
    assigned_population: str = "UNKNOWN"
    posterior_probabilities: Dict[str, float] = {}
    confidence: str = "UNKNOWN"
    quality_metrics: Dict[str, Any] = {}
    is_valid_for_scoring: bool = False
    model_hash: str = ""
    frozen_date: str = ""


# ─────────────────────── calibration_validation.json ─────────────────────
# Producer: scripts/benchmarking/27_real_world_calibration.py

class CalibrationValidationEntry(_LooseModel):
    trait: str
    population: str = "EUR"
    calibration_slope: float = 1.0
    intercept_deviation: float = 0.0
    r_squared: float = 1.0
    tail_5_accuracy: float = 0.0
    tail_95_accuracy: float = 0.0
    mean_absolute_error: float = 0.0
    is_well_calibrated: bool = True
    n_samples: int = 0
    low_confidence: bool = False


class CalibrationValidationReport(_LooseModel):
    global_status: str = ""
    mean_slope: float = 0.0
    mean_r2: float = 0.0
    well_calibrated: int = 0
    poorly_calibrated: int = 0
    low_confidence: int = 0
    generated_date: str = ""
    tolerances: Dict[str, float] = {}
    validations: List[CalibrationValidationEntry] = []


# ─────────────────────────── uncertainty_report.json ─────────────────────
# Producer: scripts/validation/14_uncertainty_propagation.py
#
# Note: the in-memory UncertaintyReport dataclass also carries a
# gwas_evidence_summary field, but the producer's actual json.dump() for
# uncertainty_report.json never writes it - so it's deliberately absent
# here too. Schema-ing the dataclass instead of the real emitted JSON
# would reject today's real output, exactly the class of bug 2.2 exists
# to catch (see tests/test_comprehensive_report_render.py's drifted
# fixture, fixed alongside this file).

class UncertaintyDecomposition(_LooseModel):
    total_variance: float = 0.0
    genotype_variance: float = 0.0
    ancestry_variance: float = 0.0
    effect_variance: float = 0.0
    genotype_fraction: float = 0.0
    ancestry_fraction: float = 0.0
    effect_fraction: float = 0.0


class UncertaintyResultEntry(_LooseModel):
    individual_id: str
    trait: str
    prs_point_estimate: float = 0.0
    prs_std_error: float = 0.0
    confidence_interval_95: Tuple[float, float] = (0.0, 0.0)
    confidence_interval_68: Tuple[float, float] = (0.0, 0.0)
    uncertainty_score: float = 0.0
    decomposition: UncertaintyDecomposition = UncertaintyDecomposition()
    n_snps_with_genotype: int = 0
    n_snps_with_effect_se: int = 0


class UncertaintyReport(_LooseModel):
    global_uncertainty_score: float = 0.0
    method: str = "three_layer_variance_propagation"
    genotype_quality_summary: Dict[str, Any] = {}
    ancestry_entropy: float = 0.0
    results: List[UncertaintyResultEntry] = []


# ───────────────────────────────  Helpers ─────────────────────────────────

def validate_and_write_json(model_instance: BaseModel, path: Union[str, Path]) -> None:
    """Write a pydantic model instance as JSON.

    The model is already valid by construction (pydantic validates at
    instantiation time) - this just standardizes the write side (matching
    the indent=2 formatting the hand-built dict producers already used)
    across all producers.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(model_instance.model_dump(mode="json"), fh, indent=2)


def load_and_validate_json(model_cls: Type[ModelT], path: Union[str, Path]) -> ModelT:
    """Load a JSON file and validate it against a pydantic model.

    Raises pydantic.ValidationError identifying exactly which field/type
    is wrong, instead of a downstream `.get(key, default)` silently
    degrading to a default value.
    """
    path = Path(path)
    with open(path) as fh:
        data = json.load(fh)
    return model_cls.model_validate(data)
