"""
Jinja2 environment + top-level document orchestration for comprehensive_report.py
(IMPROVEMENT_PLAN.md 1.6, Phase 2).

Phase 2 migrates only the outer document *shell* (head/style/header/footer/
script) to a template - every section builder still returns a pre-built HTML
string, passed into the template as one opaque `sections_html` variable, to
keep this phase's risk isolated to ~30 lines of shell markup. Section-level
migration happens in Phase 3.
"""

from datetime import datetime
from pathlib import Path
from typing import Dict

from jinja2 import Environment, FileSystemLoader, select_autoescape

from utils.constants import PIPELINE_VERSION

TEMPLATES_DIR = Path(__file__).resolve().parents[3] / "templates"
STATIC_DIR = TEMPLATES_DIR / "static"

_env = Environment(
    loader=FileSystemLoader(str(TEMPLATES_DIR)),
    autoescape=select_autoescape(enabled_extensions=()),  # sections are pre-rendered HTML, not escaped
    trim_blocks=True,
    lstrip_blocks=True,
)


def _read_static(name: str) -> str:
    return (STATIC_DIR / name).read_text(encoding="utf-8")


def render_partial(name: str, **context) -> str:
    """Render one section partial from templates/partials/{name}.

    Used by every build_*_section function in sections.py/comprehensive_report.py
    (IMPROVEMENT_PLAN.md 1.6 Phase 3) instead of hand-built f-strings.
    """
    return _env.get_template(f"partials/{name}").render(**context)


def build_html_report(lang: str, data: Dict, sample_id: str, sections_html: str,
                       reference_coverage_banner_html: str, ui: dict) -> str:
    """Render the full HTML document shell around already-built sections_html."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M UTC")
    pop = data["ancestry"].get("assigned_population", "EUR")
    integrity_score = data["integrity"].get("scientific_integrity_score", 0)

    template = _env.get_template("base.html.j2")
    return template.render(
        lang=lang,
        ui=ui,
        sample_id=sample_id,
        pop=pop,
        integrity_score_fmt=f"{integrity_score:.0f}",
        now=now,
        pipeline_version=PIPELINE_VERSION,
        css=_read_static("report.css"),
        js=_read_static("report.js"),
        reference_coverage_banner_html=reference_coverage_banner_html,
        sections_html=sections_html,
    )
