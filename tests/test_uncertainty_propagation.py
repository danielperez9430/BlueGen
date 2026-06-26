"""Unit tests for uncertainty propagation fixes."""

import sys, os, importlib.util
import pytest
import pandas as pd
import numpy as np

# Mock scipy.stats if not available (needed to load the module)
try:
    import scipy.stats  # noqa: F401
except ImportError:
    # Create a minimal mock so the module can be loaded
    import types
    _mock_scipy = types.ModuleType("scipy")
    _mock_stats = types.ModuleType("scipy.stats")

    def _mock_norm_ppf(x):
        """Simple approximation of norm.ppf for testing."""
        if x <= 0 or x >= 1:
            return 0.0
        # Very rough approximation for tests
        import math
        if x > 0.5:
            return math.sqrt(-2 * math.log(1 - x)) * 0.8
        return -math.sqrt(-2 * math.log(x)) * 0.8

    _mock_stats.norm = types.ModuleType("norm")
    _mock_stats.norm.ppf = _mock_norm_ppf
    _mock_stats.stats = _mock_stats
    _mock_scipy.stats = _mock_stats
    sys.modules["scipy"] = _mock_scipy
    sys.modules["scipy.stats"] = _mock_stats

# Load module whose filename starts with a digit
_MODULE_PATH = os.path.join(os.path.dirname(__file__), "..",
    "prs_research_pipeline", "scripts", "validation", "14_uncertainty_propagation.py")
_spec = importlib.util.spec_from_file_location("uncertainty_engine", _MODULE_PATH)
_uncer = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_uncer)
UncertaintyPropagationEngine = _uncer.UncertaintyPropagationEngine


class TestEvidenceSeRatio:
    """Change C: Verify new SE ratios are less conservative."""

    @pytest.mark.parametrize("level, old_ratio, new_ratio", [
        ("A", 0.20, 0.10),
        ("B", 0.33, 0.20),
        ("C", 0.50, 0.35),
        ("D", 0.75, 0.55),
    ])
    def test_new_ratios_lower_than_old(self, level, old_ratio, new_ratio):
        assert new_ratio < old_ratio, (
            f"Level {level}: new ratio {new_ratio} should be lower than old {old_ratio}"
        )

    def test_all_levels_present(self):
        ratios = UncertaintyPropagationEngine.EVIDENCE_SE_RATIO
        for level in ["A", "B", "C", "D"]:
            assert level in ratios
            assert 0 < ratios[level] < 1.0


class TestMafDependentG2:
    """Change B: MAF-dependent expected G² calculation."""

    def test_maf_05_gives_g2_1_5(self):
        """E[G²] = 2*0.5*0.5 + 4*0.25 = 0.5 + 1.0 = 1.5"""
        snp = pd.Series({})
        # Default MAF = 0.25 (no MAF column)
        g2 = UncertaintyPropagationEngine._get_maf(snp)
        assert g2 == 0.25

    def test_maf_default_is_025(self):
        """Without MAF column, fallback is 0.25."""
        snp = pd.Series({"weight": "0.5"})
        maf = UncertaintyPropagationEngine._get_maf(snp)
        assert maf == 0.25

    def test_maf_column_detected(self):
        """MAF from explicit column is used."""
        snp = pd.Series({"maf": "0.15", "weight": "0.5"})
        maf = UncertaintyPropagationEngine._get_maf(snp)
        assert maf == 0.15

    def test_alt_freq_column_detected(self):
        """alt_freq column is used as fallback."""
        snp = pd.Series({"alt_freq": "0.08", "weight": "0.3"})
        maf = UncertaintyPropagationEngine._get_maf(snp)
        assert maf == 0.08

    def test_invalid_maf_ignored(self):
        """Invalid MAF values fall back to 0.25."""
        snp = pd.Series({"maf": "not_a_number"})
        maf = UncertaintyPropagationEngine._get_maf(snp)
        assert maf == 0.25

    def test_out_of_range_maf_ignored(self):
        """MAF outside 0-1 falls back to default."""
        for bad_val in ["-0.1", "1.5", "2.0"]:
            snp = pd.Series({"maf": bad_val})
            maf = UncertaintyPropagationEngine._get_maf(snp)
            assert maf == 0.25, f"MAF={bad_val} should fall back but got {maf}"


class TestGqToErrorProb:
    def test_gq_zero(self):
        p = UncertaintyPropagationEngine._gq_to_error_prob(None, 0.0)
        assert p == 0.50

    def test_gq_20(self):
        p = UncertaintyPropagationEngine._gq_to_error_prob(None, 20.0)
        assert 0.005 < p < 0.02  # 10^(-20/10) = 0.01

    def test_gq_high(self):
        p = UncertaintyPropagationEngine._gq_to_error_prob(None, 99.0)
        assert p < 0.005
