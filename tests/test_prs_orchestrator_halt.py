"""Regression tests for prs.py's "required stage" halting (IMPROVEMENT_PLAN.md 4.3).

Before this fix, run_script(..., required=True) and require_output() both
printed an error but returned a value that every caller in cmd_run()
discarded - a failed critical stage (VCF->PLINK, QC, LD-pruning) let the
pipeline keep running against missing/stale files instead of stopping.
These tests pin the fixed behavior: both now halt via sys.exit() instead
of returning a value nobody checks.

Importing prs.py is safe under pytest: its venv re-exec guard only fires
when the running interpreter's sys.path lacks a "site-packages" entry
containing "venv" - false when running via ./venv/bin/python (the normal
way this test suite is run), so no os.execv happens on import.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import prs  # noqa: E402


class _FakeCompletedProcess:
    def __init__(self, returncode):
        self.returncode = returncode
        self.stdout = ""
        self.stderr = ""


def test_run_script_required_true_exits_when_script_missing():
    """ROUTES.get(key, key) falls back to the key itself for an unknown
    key, which never resolves to a real file under SCRIPTS."""
    with pytest.raises(SystemExit):
        prs.run_script("__definitely_not_a_real_routes_key__", required=True)


def test_run_script_required_false_does_not_exit_when_script_missing():
    """Same missing-script path, but required=False must still just skip."""
    rc = prs.run_script("__definitely_not_a_real_routes_key__", required=False)
    assert rc == 0


def test_run_script_required_true_exits_on_nonzero_exit(monkeypatch):
    monkeypatch.setattr(prs.subprocess, "run", lambda *a, **k: _FakeCompletedProcess(1))
    with pytest.raises(SystemExit):
        prs.run_script("prs_plink_score", required=True)


def test_run_script_required_true_does_not_exit_on_success(monkeypatch):
    monkeypatch.setattr(prs.subprocess, "run", lambda *a, **k: _FakeCompletedProcess(0))
    rc = prs.run_script("prs_plink_score", required=True)
    assert rc == 0


def test_run_script_required_true_does_not_exit_on_leakage_gate(monkeypatch):
    """Exit code 77 (leakage gate blocked) is documented as non-fatal even
    when required=True - must not be treated as a halting failure."""
    monkeypatch.setattr(prs.subprocess, "run", lambda *a, **k: _FakeCompletedProcess(77))
    rc = prs.run_script("prs_plink_score", required=True)
    assert rc == 77


def test_require_output_exits_when_missing():
    with pytest.raises(SystemExit):
        prs.require_output("definitely/not/a/real/path.bed", "label", "stage_x")


def test_require_output_returns_true_when_present(monkeypatch):
    monkeypatch.setattr(prs, "exists", lambda p: True)
    assert prs.require_output("whatever", "label") is True
