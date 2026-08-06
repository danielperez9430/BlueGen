"""Unit tests for bluegen.schemas (IMPROVEMENT_PLAN.md 2.2).

Each model gets: a valid-fixture-parses test, a couple of invalid-fixture-
rejects tests (missing required field, wrong type), and a regression test
using a small hand-copied fixture matching the real producer's actual
current output shape (not the user's real genome-derived data) - the
check that would have caught tests/test_comprehensive_report_render.py's
drifted field names immediately instead of silently degrading via
.get(key, default).
"""

import sys
from pathlib import Path

import pytest
from pydantic import ValidationError

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline"))

from bluegen.schemas import (
    PRSResult, PRSResultEntry,
    AncestryModel,
    CalibrationValidationReport, CalibrationValidationEntry,
    UncertaintyReport, UncertaintyResultEntry, UncertaintyDecomposition,
    validate_and_write_json, load_and_validate_json,
)


# ─────────────────────────── PRS_RESULT.json ────────────────────────────

class TestPRSResult:
    def test_valid_minimal(self):
        r = PRSResult(sample_id="SAMPLE_001", prs_entries=[PRSResultEntry(trait="LDL cholesterol")])
        assert r.prs_entries[0].risk_category == "medium"

    def test_missing_trait_rejected(self):
        with pytest.raises(ValidationError):
            PRSResultEntry()

    def test_missing_sample_id_rejected(self):
        with pytest.raises(ValidationError):
            PRSResult(prs_entries=[])

    def test_wrong_type_rejected(self):
        with pytest.raises(ValidationError):
            PRSResultEntry(trait="LDL cholesterol", raw_score="not_a_number")

    def test_matches_real_producer_shape(self):
        """Fixture shaped like scripts/sss/37_prs_result_unified.py's actual output."""
        fixture = {
            "sample_id": "SAMPLE_001",
            "pipeline_version": "1.2.0",
            "generated_date": "2026-08-07 12:00 UTC",
            "result_hash": "abc123def4567890",
            "ancestry": {
                "assigned_population": "EUR",
                "confidence": "HIGH",
                "probabilities": {"EUR": 0.95, "AFR": 0.02, "EAS": 0.01, "SAS": 0.01, "AMR": 0.01},
            },
            "metadata": {
                "n_traits": 1,
                "n_sources_available": 5,
                "computation_method": "PLINK --score (dosage-weighted)",
                "prs_formula": "PRS = Σ(βⱼ × Gᵢⱼ)",
                "pipeline_version": "1.2.0",
                "consolidation_note": "Single Source of Scientific Truth",
                "reference_coverage": "genome_wide",
            },
            "prs_entries": [{
                "trait": "LDL cholesterol", "raw_score": 1.234, "pca_adjusted_score": 1.1,
                "ancestry_adjusted_score": 1.1, "population_percentile": 82.0,
                "population_zscore": 0.92, "uncertainty_score": 0.3,
                "ci_95_lower": 0.5, "ci_95_upper": 1.7, "risk_category": "high",
                "assigned_population": "EUR", "calibration_mu": 0.0, "calibration_sigma": 1.0,
                "n_snps_used": 8, "n_snps_total": 10, "computation_method": "PLINK --score",
            }],
        }
        r = PRSResult.model_validate(fixture)
        assert r.prs_entries[0].population_zscore == 0.92
        assert r.prs_entries[0].population_percentile == 82.0
        # the field names this schema exists to enforce - not z_score/percentile/evidence_level
        assert "population_zscore" in fixture["prs_entries"][0]
        assert "z_score" not in fixture["prs_entries"][0]


# ─────────────────────────── ANCESTRY_MODEL.json ────────────────────────

class TestAncestryModel:
    def test_valid_defaults(self):
        m = AncestryModel()
        assert m.assigned_population == "UNKNOWN"

    def test_wrong_type_rejected(self):
        with pytest.raises(ValidationError):
            AncestryModel(n_pcs="twenty")

    def test_matches_producer_a_shape(self):
        """scripts/sss/39_ancestry_model_unified.py's quality_metrics shape."""
        fixture = {
            "method": "PCA_ENSEMBLE_V2", "reference_panel": "1000 Genomes Phase 3 (all autosomes)",
            "n_pcs": 20, "n_reference_samples": 2504,
            "super_populations": ["EUR", "AFR", "EAS", "SAS", "AMR"],
            "assigned_population": "EUR",
            "posterior_probabilities": {"EUR": 0.9, "AFR": 0.05, "EAS": 0.02, "SAS": 0.02, "AMR": 0.01},
            "confidence": "HIGH",
            "quality_metrics": {
                "max_probability": 0.9, "entropy": 0.4, "distance_ratio": 3.2, "n_pcs_used": 20,
            },
            "is_valid_for_scoring": True,
            "model_hash": "deadbeef01234567",
            "frozen_date": "2026-08-07 12:00 UTC",
        }
        m = AncestryModel.model_validate(fixture)
        assert m.quality_metrics["max_probability"] == 0.9

    def test_matches_producer_b_shape(self):
        """scripts/stages/pca_ancestry_classifier.py's differently-shaped quality_metrics."""
        fixture = {
            "method": "PCA_ENSEMBLE_V2", "reference_panel": "1000 Genomes Phase 3 (all autosomes)",
            "n_pcs": 20, "n_reference_samples": 2504,
            "super_populations": ["EUR", "AFR", "EAS", "SAS", "AMR"],
            "assigned_population": "EUR",
            "posterior_probabilities": {"EUR": 0.9, "AFR": 0.05, "EAS": 0.02, "SAS": 0.02, "AMR": 0.01},
            "confidence": "HIGH",
            "quality_metrics": {"confidence_ratio": 6.5, "knn_mean_distance": 1.23},
            "is_valid_for_scoring": True,
            "model_hash": "deadbeef01234567",
            "frozen_date": "2026-08-07 12:00 UTC",
        }
        m = AncestryModel.model_validate(fixture)
        assert m.quality_metrics["confidence_ratio"] == 6.5


# ─────────────────────── calibration_validation.json ────────────────────

class TestCalibrationValidationReport:
    def test_valid_minimal(self):
        r = CalibrationValidationReport(validations=[CalibrationValidationEntry(trait="LDL cholesterol")])
        assert r.validations[0].population == "EUR"

    def test_missing_trait_rejected(self):
        with pytest.raises(ValidationError):
            CalibrationValidationEntry()

    def test_wrong_type_rejected(self):
        with pytest.raises(ValidationError):
            CalibrationValidationEntry(trait="LDL cholesterol", is_well_calibrated=["not", "a", "bool"])

    def test_matches_real_producer_shape(self):
        """Fixture shaped like scripts/benchmarking/27_real_world_calibration.py's output."""
        fixture = {
            "global_status": "GOOD", "mean_slope": 1.02, "mean_r2": 0.91,
            "well_calibrated": 8, "poorly_calibrated": 2, "low_confidence": 1,
            "generated_date": "2026-08-07 12:00 UTC",
            "tolerances": {"slope": 0.15, "r2": 0.80},
            "validations": [{
                "trait": "LDL cholesterol", "population": "EUR",
                "calibration_slope": 1.02, "intercept_deviation": 0.02, "r_squared": 0.98,
                "tail_5_accuracy": 1.0, "tail_95_accuracy": 1.0, "mean_absolute_error": 0.01,
                "is_well_calibrated": True, "n_samples": 0, "low_confidence": False,
            }],
        }
        r = CalibrationValidationReport.model_validate(fixture)
        assert r.validations[0].is_well_calibrated is True


# ─────────────────────────── uncertainty_report.json ─────────────────────

class TestUncertaintyReport:
    def test_valid_minimal(self):
        r = UncertaintyReport(results=[
            UncertaintyResultEntry(individual_id="SAMPLE_001", trait="LDL cholesterol")
        ])
        assert r.results[0].decomposition.total_variance == 0.0

    def test_missing_trait_rejected(self):
        with pytest.raises(ValidationError):
            UncertaintyResultEntry(individual_id="SAMPLE_001")

    def test_matches_real_producer_shape(self):
        """Fixture shaped like scripts/validation/14_uncertainty_propagation.py's output.

        Note: no gwas_evidence_summary key - the real producer's json.dump()
        never writes that field despite it existing on the in-memory
        dataclass, so the schema deliberately doesn't require it either.
        """
        fixture = {
            "global_uncertainty_score": 0.42, "method": "three_layer_variance_propagation",
            "genotype_quality_summary": {"mean_uncertainty": 0.01, "n_snps": 0},
            "ancestry_entropy": 0.35,
            "results": [{
                "individual_id": "SAMPLE_001", "trait": "LDL cholesterol",
                "prs_point_estimate": 1.234, "prs_std_error": 0.135,
                "confidence_interval_95": [0.97, 1.5], "confidence_interval_68": [1.1, 1.37],
                "uncertainty_score": 1.0,
                "decomposition": {
                    "total_variance": 0.0182, "genotype_variance": 0.0003,
                    "ancestry_variance": 0.0, "effect_variance": 0.0179,
                    "genotype_fraction": 0.016, "ancestry_fraction": 0.0, "effect_fraction": 0.98,
                },
                "n_snps_with_genotype": 8, "n_snps_with_effect_se": 10,
            }],
        }
        r = UncertaintyReport.model_validate(fixture)
        assert r.results[0].decomposition.effect_fraction == 0.98

    def test_genotype_variance_string_quirk_tolerated_on_read(self):
        """Real producer bug: numpy float32 -> json.dump(default=str) silently
        stringifies decomposition.genotype_variance in some rows. Reading
        should coerce it back to a float rather than reject real output
        (the write-side fix belongs to the producer wiring in 2.2 Step 4)."""
        entry = UncertaintyResultEntry.model_validate({
            "individual_id": "SAMPLE_001", "trait": "LDL cholesterol",
            "decomposition": {"genotype_variance": "0.000296"},
        })
        assert entry.decomposition.genotype_variance == pytest.approx(0.000296)
        assert isinstance(entry.decomposition.genotype_variance, float)


# ───────────────────────────────  Helpers ─────────────────────────────────

class TestHelpers:
    def test_validate_and_write_json_roundtrip(self, tmp_path):
        r = PRSResult(sample_id="SAMPLE_001", prs_entries=[PRSResultEntry(trait="LDL cholesterol")])
        out = tmp_path / "PRS_RESULT.json"
        validate_and_write_json(r, out)
        loaded = load_and_validate_json(PRSResult, out)
        assert loaded.sample_id == "SAMPLE_001"
        assert loaded.prs_entries[0].trait == "LDL cholesterol"

    def test_load_and_validate_json_raises_on_missing_required_field(self, tmp_path):
        import json
        out = tmp_path / "bad.json"
        out.write_text(json.dumps({"sample_id": "SAMPLE_001", "prs_entries": [{"raw_score": 1.0}]}))
        with pytest.raises(ValidationError):
            load_and_validate_json(PRSResult, out)

    def test_extra_fields_preserved_not_rejected(self):
        r = PRSResult.model_validate({
            "sample_id": "SAMPLE_001", "prs_entries": [], "some_future_field": 42,
        })
        assert r.model_dump()["some_future_field"] == 42
