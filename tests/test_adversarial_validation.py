"""Unit tests for the adversarial stress-test aggregation
(scripts/publication/43_adversarial_prs_validation.py).

Regression coverage for a bug where `critical_findings` collected every
CRITICAL-severity test regardless of whether it actually failed
(`critical = [r.test_id for r in results if r.severity == "CRITICAL"]`),
instead of only the ones that failed. Effect: the HTML report's red
"Critical Findings" banner (comprehensive_report.py) listed tests that had
passed (is_robust=True) right next to a table showing them as "Robust" -
a visible self-contradiction in the same report section.
"""

import importlib.util
import os

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
