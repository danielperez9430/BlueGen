"""Unit tests for the publication evidence pack consistency checks
(scripts/publication/45_publication_evidence_pack.py).

Regression coverage for the "Critical failures" dimension: it used to read
failure_map.n_critical, the size of a static, hardcoded catalog of known
failure modes (44_failure_mode_map.py's FAILURE_MODES list) that is always
7 regardless of pipeline state - comparing it against "expected: 0" could
never pass by construction, and it silently masked the real (now-fixed)
signal of how many CRITICAL-severity adversarial tests are actually
failing right now.
"""

import json
import importlib.util
import os

_MODULE_PATH = os.path.join(os.path.dirname(__file__), "..",
    "prs_research_pipeline", "scripts", "publication", "45_publication_evidence_pack.py")
_spec = importlib.util.spec_from_file_location("evidence_pack_module", _MODULE_PATH)
_pack = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_pack)
PublicationEvidencePack = _pack.PublicationEvidencePack


def _write_json(path, data):
    with open(path, "w") as fh:
        json.dump(data, fh)


def _make_inputs(tmp_path, critical_findings, n_critical_catalog=7):
    """Writes minimal, otherwise-passing input JSONs and returns their paths."""
    prs_core = tmp_path / "prs_core_definition.json"
    _write_json(prs_core, {"formula": "PRS = Σ(βⱼ × Gᵢⱼ)", "definition_hash": "abc123"})

    ancestry = tmp_path / "ANCESTRY_MODEL.json"
    _write_json(ancestry, {"method": "PCA_ENSEMBLE_V2", "is_valid_for_scoring": True})

    integrity = tmp_path / "FINAL_SCIENTIFIC_SCORE.json"
    _write_json(integrity, {"scientific_integrity_score": 84.4})

    adversarial = tmp_path / "adversarial_validation_report.json"
    _write_json(adversarial, {"overall_robustness_score": 85.0,
                               "critical_findings": critical_findings})

    # Static catalog count deliberately different from critical_findings,
    # to prove the check reads from adversarial, not from this file.
    failure_map = tmp_path / "failure_mode_map.json"
    _write_json(failure_map, {"n_critical": n_critical_catalog})

    return dict(prs_core_json=str(prs_core), ancestry_json=str(ancestry),
                integrity_json=str(integrity), adversarial_json=str(adversarial),
                failure_map_json=str(failure_map))


def test_critical_failures_dimension_reads_active_adversarial_findings(tmp_path):
    pack = PublicationEvidencePack(output_dir=str(tmp_path / "out"))
    inputs = _make_inputs(tmp_path, critical_findings=["LD_DISRUPT_SEVERE"])

    result = pack.generate(**inputs)

    crit_check = next(c for c in result.consistency if c.dimension == "Critical failures")
    assert crit_check.observed == "1"
    assert crit_check.match is False
    # Must not be the static catalog's count (7), which was the old bug.
    assert crit_check.observed != "7"


def test_critical_failures_dimension_passes_when_no_active_findings(tmp_path):
    pack = PublicationEvidencePack(output_dir=str(tmp_path / "out"))
    # Static catalog still has 7 "known" critical failure modes documented -
    # that must not block this check once nothing is actively failing.
    inputs = _make_inputs(tmp_path, critical_findings=[], n_critical_catalog=7)

    result = pack.generate(**inputs)

    crit_check = next(c for c in result.consistency if c.dimension == "Critical failures")
    assert crit_check.observed == "0"
    assert crit_check.match is True
