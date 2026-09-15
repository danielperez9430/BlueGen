"""End-to-end render regression test for comprehensive_report.py
(IMPROVEMENT_PLAN.md 1.6 safety net).

Before this test, comprehensive_report.py had 90 unit tests on isolated
helper functions but nothing exercising the full build_html_report() path -
nothing would catch a broken section wire-up, a missing section, or (once
the Jinja2 migration lands) an unrendered template leftover. This test uses
only small hand-written synthetic fixtures (safe to commit - no real genome
data) and must keep passing, unchanged, through every phase of that refactor.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__),
    "..", "prs_research_pipeline", "scripts", "publication"
))

from comprehensive_report import build_html_report  # noqa: E402


def _fake_data():
    # Field names here must match bluegen.schemas.PRSResultEntry / the real
    # producer (scripts/sss/37_prs_result_unified.py) - population_zscore/
    # population_percentile, no evidence_level - not the PGS-section's
    # unrelated z_score/percentile/evidence_level fields (IMPROVEMENT_PLAN.md 2.2).
    prs_entries = [
        {
            "trait": "Lactose intolerance", "population_zscore": 2.3, "risk_category": "high",
            "population_percentile": 95, "n_snps_used": 3, "n_snps_total": 3,
            "uncertainty_score": 0.2,
        },
        {
            "trait": "Morning chronotype (early bird)", "population_zscore": 1.8, "risk_category": "high",
            "population_percentile": 90, "n_snps_used": 2, "n_snps_total": 2,
            "uncertainty_score": 0.3,
        },
        {
            "trait": "Vitamin D metabolism", "population_zscore": -0.5, "risk_category": "low",
            "population_percentile": 30, "n_snps_used": 8, "n_snps_total": 12,
            "uncertainty_score": 0.4,
        },
    ]
    return {
        "prs_result": {
            "sample_id": "TEST_SAMPLE",
            "prs_entries": prs_entries,
            "metadata": {"reference_coverage": "genome_wide"},
        },
        "ancestry": {
            "assigned_population": "EUR",
            "posterior_probabilities": {"EUR": 0.98, "AFR": 0.01, "EAS": 0.01},
            "confidence": "HIGH",
            "n_reference_samples": 2504,
        },
        "validation": {"checks": []},
        "integrity": {"scientific_integrity_score": 85.0, "integrity_category": "RESEARCH_GRADE"},
        "benchmark": {},
        "adversarial": {},
        "failure_map": {},
        "leakage": {},
        "quality_delta": {},
        "uncertainty_report": {"results": []},
        "calibration_report": {},
        "calibration_validation": {"validations": []},
        "gwas_consortium": {},
        "portability": {},
        "reproducibility": {},
        "consistency": {},
        "leakage_audit": {},
        "snp_universe": {},
        "pgs_calibration": {"all_entries": []},
        "clinvar": {},
        "pharmgkb": {},
        "deep_ancestry": {},
        "_cal_lookup": {}, "_uncert_lookup": {}, "_evidence_lookup": {},
        "_recommendation_lookup": {}, "_polarity_inverted": set(),
        "_pgs_coverage_lookup": {},
    }


UNCONDITIONAL_SECTION_IDS = [
    "top_findings", "summary", "ancestry", "prs", "uncertainty_decomp",
    "variants", "calibration", "clinvar", "clinical_actionability",
    "portability", "validation", "gwas_consortium", "benchmark",
    "adversarial", "failure_map", "leakage", "consistency", "integrity",
    "radar", "reproducibility", "methodology", "limitations",
]


class TestBuildHtmlReportRendersWithoutTemplateLeftovers:
    @pytest.mark.parametrize("lang", ["en", "es"])
    def test_no_unrendered_jinja_or_fstring_artifacts(self, lang):
        html = build_html_report(lang, _fake_data(), "TEST_SAMPLE")
        assert "{{" not in html
        assert "}}" not in html
        assert "{%" not in html

    @pytest.mark.parametrize("lang", ["en", "es"])
    def test_is_a_well_formed_document(self, lang):
        html = build_html_report(lang, _fake_data(), "TEST_SAMPLE")
        assert html.startswith("<!DOCTYPE html>")
        assert "</html>" in html
        assert f'<html lang="{lang}">' in html


class TestBuildHtmlReportSectionWiring:
    @pytest.mark.parametrize("lang", ["en", "es"])
    @pytest.mark.parametrize("section_id", UNCONDITIONAL_SECTION_IDS)
    def test_every_unconditional_section_id_present(self, lang, section_id):
        html = build_html_report(lang, _fake_data(), "TEST_SAMPLE")
        assert f'id="{section_id}"' in html

    def test_sample_id_appears_in_header(self):
        html = build_html_report("en", _fake_data(), "TEST_SAMPLE")
        assert "TEST_SAMPLE" in html


class TestEvidenceLevelLegend:
    """The A-D evidence letters must be explained where they first appear
    (Top Findings) and in Variant-Level Detail, in both languages, with the
    same definitions that head data/snp_database.csv (user request 2026-09-15)."""

    @pytest.mark.parametrize("lang", ["en", "es"])
    def test_legend_present_in_both_sections(self, lang):
        html = build_html_report(lang, _fake_data(), "TEST_SAMPLE")
        assert html.count('class="evidence-legend"') >= 2
        title = "Niveles de evidencia" if lang == "es" else "Evidence levels"
        assert title in html
        for letter in "ABCD":
            assert f'font-weight:700">{letter}</span>' in html

    def test_definitions_match_the_panel_preamble(self):
        import re
        from pathlib import Path
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "prs_research_pipeline", "scripts", "publication"))
        from report.interpretations import EVIDENCE_LEVELS
        preamble = Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "data" / "snp_database.csv"
        head = "".join(line for line in preamble.read_text().splitlines(keepends=True)[:6] if line.startswith("#"))
        panel = dict(re.findall(r"([ABCD])=([^,\n]+)", head))
        assert set(panel) == {"A", "B", "C", "D"}
        legend = dict(EVIDENCE_LEVELS["en"])
        # same spirit, checked by the key term of each definition
        assert "5×10" in legend["A"] and "GWAS" in panel["A"]
        assert "candidate" in legend["B"].lower() and "candidate" in panel["B"].lower()
        assert "single study" in legend["C"].lower() and "single study" in panel["C"].lower()
        assert "mechanistic" in legend["D"].lower() and "mechanistic" in panel["D"].lower()
