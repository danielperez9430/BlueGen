"""Guards against the pipeline version drifting across files again.

utils.constants.PIPELINE_VERSION is the single source of truth (see TIER 0.1
in IMPROVEMENT_PLAN.md). Before this test existed, README.md, config.yaml,
prs.py, and six different scripts each hardcoded their own independent
version string - resulting in v1.0.0, v1.1.0, v2.0.0, v6.0.0, v7.0.0 and
v9.0.0 all appearing in the same generated report.
"""

import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "prs_research_pipeline" / "scripts"

sys.path.insert(0, str(SCRIPTS_DIR))
from utils.constants import PIPELINE_VERSION


def test_root_readme_badge_matches():
    text = (REPO_ROOT / "README.md").read_text()
    m = re.search(r"badge/version-([0-9]+\.[0-9]+\.[0-9]+)-", text)
    assert m, "version badge not found in README.md"
    assert m.group(1) == PIPELINE_VERSION


def test_pipeline_readme_table_matches():
    text = (REPO_ROOT / "prs_research_pipeline" / "README.md").read_text()
    versions_in_table = set(re.findall(r"^\|\s*([0-9]+\.[0-9]+\.[0-9]+)\s*\|", text, re.MULTILINE))
    assert versions_in_table, "no version row found in prs_research_pipeline/README.md"
    assert versions_in_table == {PIPELINE_VERSION}, (
        f"prs_research_pipeline/README.md lists version(s) {versions_in_table}, "
        f"expected only {{{PIPELINE_VERSION!r}}}"
    )


def test_config_yaml_matches():
    config = yaml.safe_load((REPO_ROOT / "prs_research_pipeline" / "config.yaml").read_text())
    assert config["pipeline"]["version"] == PIPELINE_VERSION


def test_no_other_hardcoded_pipeline_version_literals():
    """Every pipeline_version field must reference the constant, not a literal."""
    pattern = re.compile(r'pipeline_version["\']?\s*[:=]\s*["\']([0-9]+\.[0-9]+\.[0-9]+)["\']')
    offenders = []
    for py_file in SCRIPTS_DIR.rglob("*.py"):
        if py_file.name == "constants.py":
            continue
        for m in pattern.finditer(py_file.read_text()):
            offenders.append(f"{py_file.relative_to(REPO_ROOT)}: {m.group(0)}")
    assert not offenders, "hardcoded pipeline_version literal(s) found:\n" + "\n".join(offenders)


def test_prs_py_uses_the_constant_not_a_literal():
    text = (REPO_ROOT / "prs.py").read_text()
    assert "PIPELINE_VERSION" in text
    assert not re.search(r'BlueGen v[0-9]+\.[0-9]+\.[0-9]+["{]', text.replace("v{PIPELINE_VERSION}", "")), (
        "prs.py appears to hardcode a version string instead of using PIPELINE_VERSION"
    )
