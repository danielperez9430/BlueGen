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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline"))

from bluegen.calibration import PopulationCalibrationV2, PopulationDistribution


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
        distributions = {trait: {"EUR": dist}}
    else:
        # calibrate_sample() bails out early if distributions is entirely
        # empty, so seed an unrelated trait to keep it non-empty.
        distributions = {"Other trait": {"EUR": _dist("Other trait", 200)}}

    with tempfile.TemporaryDirectory() as tmpdir:
        sample_prs = _write_sample_prs(tmpdir, trait, prs_raw)
        ancestry_json = _write_ancestry(tmpdir)
        results = calibrator.calibrate_sample(
            sample_prs=sample_prs, ancestry_json=ancestry_json,
            output_dir=tmpdir, ref_dist_dir=None, distributions=distributions,
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


def test_skewed_reference_uses_empirical_percentile_and_clips_z():
    """Few-SNP traits give discrete, strongly skewed reference distributions
    (real case: 3-SNP 'Hair color (black)' in EUR — mean 0.007, sd 0.019,
    skew 2.7 — where a carrier of every effect allele got z = +15.5). With
    |skewness| > 2 the percentile must come from the stored quantiles, z is
    clipped to ±6 and low_confidence is set (RELEASE_PLAN 3.0.3)."""
    dist = PopulationDistribution(
        trait="Hair color (black)", population="EUR", n_samples=503,
        mean=0.007, std=0.0189, median=0.0, iqr=0.0,
        percentile_5=0.0, percentile_25=0.0, percentile_75=0.0, percentile_95=0.05,
        skewness=2.681, kurtosis=8.0, shapiro_p=0.0,
    )
    r = _calibrate("Hair color (black)", 0.30, dist)
    assert r.z_score_population == 6.0
    assert r.percentile_population == 97.5
    assert r.low_confidence is True
    assert r.risk_category == "high"
    # inside the quantile range: linear between p75 (0.0 → 75) and p95 (0.05 → 95)
    r2 = _calibrate("Hair color (black)", 0.025, dist)
    assert abs(r2.percentile_population - 85.0) < 0.6
    assert r2.low_confidence is True
    # a symmetric reference is untouched by the rule
    r3 = _calibrate("Normal trait", 1.0, _dist("Normal trait", 200, mean=0.0, std=1.0))
    assert r3.low_confidence is False and abs(r3.z_score_population - 1.0) < 1e-6


def test_discrete_symmetric_reference_is_also_low_confidence():
    """Real case: 3-SNP 'Skin pigmentation' in EUR has median 0.1 and IQR 0
    (most Europeans carry the same alleles) with skewness only −0.24; the
    normal z for a 0.26 scorer was +5.8. IQR == 0 triggers the same handling."""
    dist = PopulationDistribution(
        trait="Skin pigmentation", population="EUR", n_samples=503,
        mean=0.1039, std=0.0265, median=0.1, iqr=0.0,
        percentile_5=0.05, percentile_25=0.1, percentile_75=0.1, percentile_95=0.1417,
        skewness=-0.2388, kurtosis=2.26, shapiro_p=0.0,
    )
    r = _calibrate("Skin pigmentation", 0.2583, dist)
    assert r.low_confidence is True
    assert r.percentile_population == 97.5
    assert abs(r.z_score_population - 5.83) < 0.05   # under the ±6 clip, reported as is
    r2 = _calibrate("Skin pigmentation", 0.1, dist)  # the modal value
    assert r2.low_confidence is True and 25.0 <= r2.percentile_population <= 75.0
