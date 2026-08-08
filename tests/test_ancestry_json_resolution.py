"""Regression test for prs.py's ancestry-classification file resolution
(IMPROVEMENT_PLAN.md TIER 3.5).

cmd_run's population calibration used to fall back to
"pca/ancestry_inference.json" - a stale file from the pre-Phase-6
allele-frequency ancestry method that no current script writes, and whose
JSON has no top-level "assigned_population" key. bluegen/calibration.py's
calibrate_sample() does `ancestry.get("assigned_population", "EUR")`, so
reading that legacy file silently calibrated every sample as EUR regardless
of their real ancestry. resolve_ancestry_json() must only ever return a
path from the canonical producer(s), matching the same policy already
documented in scripts/sss/39_ancestry_model_unified.py.

Importing prs.py is safe under pytest: its venv re-exec guard only fires
when the running interpreter's sys.path lacks a "site-packages" entry
containing "venv" - false when running via ./venv/bin/python.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import prs  # noqa: E402


def test_prefers_current_pca_classifier_output(monkeypatch):
    present = {"pca/ancestry_classification.json", "ancestry/classification_report.json"}
    monkeypatch.setattr(prs, "exists", lambda p: p in present)
    assert prs.resolve_ancestry_json() == "pca/ancestry_classification.json"


def test_falls_back_to_documented_alternate_path(monkeypatch):
    present = {"ancestry/classification_report.json"}
    monkeypatch.setattr(prs, "exists", lambda p: p in present)
    assert prs.resolve_ancestry_json() == "ancestry/classification_report.json"


def test_never_resolves_the_stale_legacy_file(monkeypatch):
    """Even if the deprecated file is the only ancestry-shaped file on disk,
    it must not be picked - its schema silently defaults calibration to EUR."""
    present = {"pca/ancestry_inference.json"}
    monkeypatch.setattr(prs, "exists", lambda p: p in present)
    assert prs.resolve_ancestry_json() is None


def test_returns_none_when_classifier_never_ran(monkeypatch):
    monkeypatch.setattr(prs, "exists", lambda p: False)
    assert prs.resolve_ancestry_json() is None
