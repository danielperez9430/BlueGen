"""Unit tests for the adversarial stress-test aggregation
(scripts/publication/43_adversarial_prs_validation.py).

Regression coverage for two bugs:

1. `critical_findings` collected every CRITICAL-severity test regardless of
   whether it actually failed (`critical = [r.test_id for r in results if
   r.severity == "CRITICAL"]`), instead of only the ones that failed.
   Effect: the HTML report's red "Critical Findings" banner
   (comprehensive_report.py) listed tests that had passed (is_robust=True)
   right next to a table showing them as "Robust" - a visible
   self-contradiction in the same report section.

2. `_test_ld_disruption`'s VIF computation was `(variance * inflation) /
   variance`, which algebraically always equals exactly `inflation`
   (1.3/2.0/3.0) regardless of any real PRS data - not a stress test, a
   constant. Combined with `robust = vif < 2.0`, LD_DISRUPT_MODERATE and
   LD_DISRUPT_SEVERE were guaranteed to report "vulnerable" on every run,
   forever, independent of the actual pipeline's robustness. Replaced with
   a real simulation (perturb z-scores with noise sized to the target
   variance inflation, measure whether trait ranking survives) matching
   the methodology already used by the other three stress tests.
"""

import importlib.util
import os

import numpy as np

_MODULE_PATH = os.path.join(os.path.dirname(__file__), "..",
    "prs_research_pipeline", "scripts", "publication", "43_adversarial_prs_validation.py")
_spec = importlib.util.spec_from_file_location("adversarial_module", _MODULE_PATH)
_adv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_adv)
AdversarialPRSValidator = _adv.AdversarialPRSValidator


def test_critical_findings_only_includes_failed_critical_tests(tmp_path):
    validator = AdversarialPRSValidator(seed=42, output_dir=str(tmp_path))
    report = validator.run_all(prs_result_json=str(tmp_path / "missing_prs_result.json"))

    critical_ids = {r.test_id for r in report.results if r.severity == "CRITICAL"}
    failed_critical_ids = {r.test_id for r in report.results
                            if r.severity == "CRITICAL" and not r.is_robust}
    passed_critical_ids = critical_ids - failed_critical_ids

    # Sanity: this fixture (default seed, synthetic fallback data) must
    # actually exercise both a passing and a failing CRITICAL test, or the
    # test isn't proving anything about the filter.
    assert passed_critical_ids, "expected at least one CRITICAL test to pass"
    assert failed_critical_ids, "expected at least one CRITICAL test to fail"

    assert set(report.critical_findings) == failed_critical_ids
    for passed_id in passed_critical_ids:
        assert passed_id not in report.critical_findings


def test_n_vulnerable_matches_non_robust_count(tmp_path):
    validator = AdversarialPRSValidator(seed=42, output_dir=str(tmp_path))
    report = validator.run_all(prs_result_json=str(tmp_path / "missing_prs_result.json"))

    assert report.n_vulnerable == sum(1 for r in report.results if not r.is_robust)
    assert report.n_robust == sum(1 for r in report.results if r.is_robust)


def _entries(zscores, percentiles=None):
    if percentiles is None:
        percentiles = [50 + 10 * z for z in zscores]
    return [{"trait": f"t{i}", "population_zscore": z, "population_percentile": p}
            for i, (z, p) in enumerate(zip(zscores, percentiles))]


def test_ld_disruption_is_not_a_constant_function_of_inflation_alone():
    """The old bug made LD_DISRUPT_{label}'s outcome depend only on the
    hardcoded inflation constant, never on the actual PRS data - so two
    wildly different datasets always produced identical stressed/detail
    values. A real stress test must not do that."""
    validator = AdversarialPRSValidator(seed=7)

    tight = validator._test_ld_disruption(_entries([0.1, -0.05, 0.08, 0.02] * 5))
    validator2 = AdversarialPRSValidator(seed=7)
    spread = validator2._test_ld_disruption(_entries([-4.0, -1.5, 0.5, 3.0, 5.5] * 4))

    tight_by_id = {r.test_id: r for r in tight}
    spread_by_id = {r.test_id: r for r in spread}

    # Same seed, same target inflation levels, very different underlying
    # variance -> the measured/stressed values must differ. Under the old
    # formula both runs would report the exact same vif (1.3/2.0/3.0) no
    # matter what these arrays contained.
    differs = any(tight_by_id[tid].stressed != spread_by_id[tid].stressed
                  for tid in tight_by_id)
    assert differs, "LD disruption results must depend on real input data"


def test_ld_disruption_mild_less_disruptive_than_severe():
    validator = AdversarialPRSValidator(seed=7)
    results = validator._test_ld_disruption(
        _entries([-4.0, -1.5, 0.5, 3.0, 5.5, -2.2, 1.1, 4.4] * 3))
    by_id = {r.test_id: r for r in results}

    # Rank-correlation (stressed) should degrade monotonically-ish as the
    # nominal inflation grows - MILD's perturbation is smallest by
    # construction (extra_sd scales with sqrt(inflation - 1)).
    assert by_id["LD_DISRUPT_MILD"].stressed >= by_id["LD_DISRUPT_SEVERE"].stressed
