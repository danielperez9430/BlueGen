"""Tests for population_calibrate_v2.py's z-score/percentile consistency.

Regression coverage for a real bug: when a trait's 1000G reference
distribution had fewer than 10 samples, percentile_population was hardcoded
to 50.0 while z_score_population was still computed normally from real
mu/sigma - an internally inconsistent pair that silently forced
risk_category to "medium" and fed a degenerate calibration slope into
27_real_world_calibration.py (see IMPROVEMENT_PLAN.md TIER 3).
"""

import csv
import json
import sys
import tempfile
from pathlib import Path

from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from prs.population_calibrate_v2 import PopulationCalibrationV2, PopulationDistribution


def _write_sample_prs(tmpdir, trait, prs_raw):
    path = Path(tmpdir) / "sample_prs.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["trait", "prs_raw"])
        w.writerow([trait, prs_raw])
    return str(path)


def _write_ancestry(tmpdir):
    path = Path(tmpdir) / "ancestry.json"
    path.write_text(json.dumps({"assigned_population": "EUR", "posterior_probabilities": {"EUR": 1.0}}))
    return str(path)


def _dist(trait, n_samples, mean=0.0, std=1.0):
    return PopulationDistribution(
        trait=trait, population="EUR", n_samples=n_samples,
        mean=mean, std=std, median=mean, iqr=std,
        percentile_5=mean - 2 * std, percentile_25=mean - std,
        percentile_75=mean + std, percentile_95=mean + 2 * std,
        skewness=0.0, kurtosis=0.0, shapiro_p=1.0,
    )


def _calibrate(trait, prs_raw, dist):
    calibrator = PopulationCalibrationV2()
    if dist is not None:
        calibrator._reference_distributions = {trait: {"EUR": dist}}
    else:
        # calibrate_sample() bails out early if _reference_distributions is
        # entirely empty, so seed an unrelated trait to keep it non-empty.
        calibrator._reference_distributions = {"Other trait": {"EUR": _dist("Other trait", 200)}}

    with tempfile.TemporaryDirectory() as tmpdir:
        sample_prs = _write_sample_prs(tmpdir, trait, prs_raw)
        ancestry_json = _write_ancestry(tmpdir)
        results = calibrator.calibrate_sample(
            sample_prs=sample_prs, ancestry_json=ancestry_json,
            output_dir=tmpdir, ref_dist_dir=None,
        )
    assert len(results) == 1
    return results[0]


def test_low_ref_samples_z_percentile_consistent():
    result = _calibrate("Low-N trait", prs_raw=-0.5, dist=_dist("Low-N trait", n_samples=6))
    expected_pctl = round(scipy_stats.norm.cdf(result.z_score_population) * 100, 1)
    assert result.percentile_population == expected_pctl
    assert result.low_confidence is True
    assert result.n_reference_samples == 6


def test_normal_case_unchanged():
    result = _calibrate("Well-supported trait", prs_raw=-0.5, dist=_dist("Well-supported trait", n_samples=200))
    expected_pctl = round(scipy_stats.norm.cdf(result.z_score_population) * 100, 1)
    assert result.percentile_population == expected_pctl
    assert result.low_confidence is False
    assert result.n_reference_samples == 200


def test_no_distribution_neutralized():
    result = _calibrate("Missing trait", prs_raw=3.0, dist=None)
    assert result.z_score_population == 0.0
    assert result.percentile_population == 50.0
    assert result.low_confidence is True
    assert result.n_reference_samples == 0


def test_risk_category_low_confidence_no_longer_forced_medium():
    # z=-0.9 -> percentile ~18.4%, below the "low" threshold (25). Before the
    # fix this was hardcoded to percentile=50.0 -> risk_category="medium"
    # regardless of the real z-score.
    result = _calibrate("Previously broken trait", prs_raw=-0.9, dist=_dist("Previously broken trait", n_samples=6))
    assert result.percentile_population < 25
    assert result.risk_category == "low"


def test_risk_category_boundaries_still_respected():
    high = _calibrate("High trait", prs_raw=2.0, dist=_dist("High trait", n_samples=200))
    assert high.percentile_population >= 75
    assert high.risk_category == "high"

    medium = _calibrate("Medium trait", prs_raw=0.0, dist=_dist("Medium trait", n_samples=200))
    assert 25 < medium.percentile_population < 75
    assert medium.risk_category == "medium"
