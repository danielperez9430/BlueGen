"""Streamlit dashboard smoke tests (RELEASE_PLAN.md 2.1.8).

Runs dashboard.py headlessly with streamlit's AppTest and visits every page
twice: once against the real prs_research_pipeline/ outputs (whatever is
present locally / in CI) and once against an empty directory, to pin the
acceptance criterion that every page degrades to an informational message
instead of a traceback when its JSON/CSV is missing.
"""

import os
from pathlib import Path

import pytest

st = pytest.importorskip("streamlit")
pytest.importorskip("plotly")
from streamlit.testing.v1 import AppTest  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DASHBOARD = REPO_ROOT / "dashboard.py"

PAGES = [
    "📊 Overview",
    "🧬 PRS Results",
    "🥗 Recommendations",
    "🩺 PGS Catalog",
    "🔬 ClinVar Pathogenic",
    "💊 Pharmacogenomics",
    "🌍 Ancestry",
    "🦴 Archaic DNA",
    "📋 Raw Data",
]


def _run(page, pipeline_dir=None, monkeypatch=None):
    if pipeline_dir is not None:
        monkeypatch.setenv("BLUEGEN_PIPELINE_DIR", str(pipeline_dir))
    else:
        monkeypatch.delenv("BLUEGEN_PIPELINE_DIR", raising=False)
    at = AppTest.from_file(str(DASHBOARD), default_timeout=60)
    at.run()
    assert not at.exception, f"initial run raised: {at.exception}"
    at.sidebar.radio[0].set_value(page).run()
    assert not at.exception, f"page {page!r} raised: {[e.value for e in at.exception]}"
    return at


@pytest.mark.parametrize("page", PAGES)
def test_every_page_renders_against_real_outputs(page, monkeypatch):
    at = _run(page, None, monkeypatch)
    assert at.title, f"page {page!r} rendered no title"


@pytest.mark.parametrize("page", PAGES)
def test_every_page_degrades_cleanly_without_data(page, tmp_path, monkeypatch):
    """Empty pipeline dir: no JSON, no CSV, no PCA files. Must not raise."""
    at = _run(page, tmp_path, monkeypatch)
    assert at.title
    # Pages that depend on pipeline outputs must say so rather than show blanks
    if page not in ("📋 Raw Data", "🌍 Ancestry", "📊 Overview"):
        assert at.info or at.warning, f"page {page!r} showed neither st.info nor st.warning without data"


def test_sidebar_version_comes_from_the_constant(monkeypatch):
    import sys
    sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts"))
    from utils.constants import PIPELINE_VERSION
    at = _run(PAGES[0], None, monkeypatch)
    captions = " ".join(c.value for c in at.sidebar.caption)
    assert f"BlueGen v{PIPELINE_VERSION}" in captions


def test_recommendations_page_lists_curated_traits(monkeypatch):
    """With the real data dir: the curated file has 47 traits; the page must
    surface them (unfiltered) even when the local PRS run scored fewer."""
    at = _run("🥗 Recommendations", None, monkeypatch)
    # Turn off the "only scored" filter and count expanders
    boxes = [c for c in at.checkbox if "Only traits scored" in c.label]
    if boxes:
        at = boxes[0].set_value(False).run()
        assert not at.exception
    assert len(at.expander) >= 40, f"expected the curated traits as expanders, got {len(at.expander)}"


def test_pipeline_dir_override_is_honoured(tmp_path, monkeypatch):
    """The env override exists for these tests; make sure it really redirects
    every loader (a leaked absolute path would silently read real outputs)."""
    at = _run("🧬 PRS Results", tmp_path, monkeypatch)
    assert any("No PRS data" in w.value for w in at.warning)
