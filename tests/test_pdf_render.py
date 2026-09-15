"""PDF output of the comprehensive report (RELEASE_PLAN 3.0.5).

The HTML now carries a print-only cover page with a table of contents whose
page numbers WeasyPrint resolves with target-counter(); sections start on a
new page and rows/cards never split. Rendered from the same synthetic
fixture as test_comprehensive_report_render.py. Skipped when WeasyPrint (or
its system libraries) is not importable, e.g. the plain CI test job; runs
inside the Docker image.
"""

import re
import sys
from pathlib import Path  # noqa: F401 (re is used by the HTML test below)

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts" / "publication"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from comprehensive_report import build_html_report  # noqa: E402
from test_comprehensive_report_render import _fake_data  # noqa: E402

try:
    import weasyprint  # noqa: F401
    HAS_WEASYPRINT = True
except Exception:  # ImportError or missing cairo/pango
    HAS_WEASYPRINT = False


@pytest.mark.parametrize("lang", ["en", "es"])
def test_html_carries_a_print_cover_and_toc_with_anchors(lang):
    html = build_html_report(lang, _fake_data(), "TEST_SAMPLE")
    assert '<section class="print-cover">' in html
    assert ("Contenido" if lang == "es" else "Contents") in html
    hrefs = re.findall(r'<li><a href="#([A-Za-z0-9_\-]+)">', html)
    assert len(hrefs) >= 5
    for h in hrefs:
        assert f'id="{h}"' in html, f"TOC entry #{h} has no matching section id"
    assert 'data-sample="TEST_SAMPLE"' in html
    # every cover label rendered (a dict key named "pop" once resolved to dict.pop)
    assert ("<th>Población</th>" if lang == "es" else "<th>Population</th>") in html


@pytest.mark.skipif(not HAS_WEASYPRINT, reason="WeasyPrint not importable here")
def test_pdf_has_cover_toc_and_multiple_pages(tmp_path):
    from weasyprint import HTML
    html = build_html_report("en", _fake_data(), "TEST_SAMPLE")
    doc = HTML(string=html, base_url=str(tmp_path)).render()
    # one page per section (break-before: page) + the cover: well above 4
    assert len(doc.pages) >= 4, f"expected a cover + several section pages, got {len(doc.pages)}"
    pdf_path = tmp_path / "report.pdf"
    doc.write_pdf(str(pdf_path))
    data = pdf_path.read_bytes()
    assert data[:5] == b"%PDF-"
    assert pdf_path.stat().st_size > 20_000
