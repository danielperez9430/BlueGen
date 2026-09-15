"""Multi-sample comparison report (RELEASE_PLAN 3.0.4): renders from
synthetic per-sample calibrated CSVs, no real data."""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts" / "publication"))

import comparison_report as cr  # noqa: E402


@pytest.fixture
def pipeline(tmp_path):
    (tmp_path / "prs" / "pgs_scores").mkdir(parents=True)
    (tmp_path / "science").mkdir()
    prs = pd.DataFrame([
        {"individual_id": s, "trait": t, "z_score_population": z, "percentile_population": p,
         "risk_category": "medium", "low_confidence": lc}
        for s, t, z, p, lc in [
            ("A", "Lactose intolerance", 0.1, 54.0, False), ("B", "Lactose intolerance", -1.4, 8.1, False),
            ("C", "Lactose intolerance", 2.3, 98.9, True),
            ("A", "Iron levels", 1.2, 88.5, False), ("B", "Iron levels", 1.1, 86.4, False),
        ]])
    prs.to_csv(tmp_path / "prs" / "population_calibrated_v2.csv", index=False)
    pgs = pd.DataFrame([
        {"individual_id": s, "pgs_id": "PGS000020", "trait": "Type 2 diabetes (T2D)", "z_score": z,
         "percentile": p, "reliable": True}
        for s, z, p in [("A", 1.8, 96.4), ("B", -0.2, 42.1), ("C", 0.4, 65.5)]])
    pgs.to_csv(tmp_path / "prs" / "pgs_scores" / "pgs_calibrated.csv", index=False)
    (tmp_path / "science" / "ANCESTRY_MODEL.json").write_text(json.dumps({"assigned_population": "EUR"}))
    return tmp_path


def test_builds_three_columns_in_both_languages(pipeline):
    prs, pgs, anc = cr.load_tables(pipeline)
    for lang in ("en", "es"):
        html = cr.render(prs, pgs, anc, lang)
        assert html.count("<th>A</th>") == 3 and "<th>B</th>" in html and "<th>C</th>" in html
        assert "Lactose intolerance" in html and "PGS000020" in html
        assert "+2.30 (99%) ‡" in html            # low-confidence marker
        assert 'class="hi">+1.80' in html and 'class="lo">-1.40' in html
        assert "—" in html                        # C has no Iron levels row
    assert "Comparación" in cr.render(prs, pgs, anc, "es")


def test_rows_sorted_by_spread_and_run_level_ancestry_marked(pipeline):
    prs, pgs, anc = cr.load_tables(pipeline)
    ctx = cr.build_context(prs, pgs, anc, "en")
    assert [r["name"] for r in ctx["prs_rows"]] == ["Lactose intolerance", "Iron levels"]
    assert ctx["anc_cells"] == ["EUR", "EUR", "EUR"] and ctx["anc_note"] == "run-level (not per sample)"


def test_cli_writes_files_and_skips_single_sample(pipeline, tmp_path, capsys):
    out = tmp_path / "reports"
    assert cr.main(["--pipeline-dir", str(pipeline), "--output-dir", str(out), "--lang", "both"]) == 0
    assert (out / "comparison_report_en.html").exists() and (out / "comparison_report_es.html").exists()
    # single sample → nothing to compare
    single = pd.read_csv(pipeline / "prs" / "population_calibrated_v2.csv")
    single[single.individual_id == "A"].to_csv(pipeline / "prs" / "population_calibrated_v2.csv", index=False)
    pgs = pd.read_csv(pipeline / "prs" / "pgs_scores" / "pgs_calibrated.csv")
    pgs[pgs.individual_id == "A"].to_csv(pipeline / "prs" / "pgs_scores" / "pgs_calibrated.csv", index=False)
    out2 = tmp_path / "reports2"
    assert cr.main(["--pipeline-dir", str(pipeline), "--output-dir", str(out2)]) == 0
    assert not out2.exists()
    assert "nothing to compare" in capsys.readouterr().out
