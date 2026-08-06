"""Unit tests for comprehensive report confidence helper functions."""

import sys
import os
import pytest

# Add the scripts/publication directory to path
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__),
    "..", "prs_research_pipeline", "scripts", "publication"
))

from comprehensive_report import (
    compute_per_trait_confidence,
    confidence_stars,
    calibration_flag,
    trust_tier,
    trust_badge,
    mini_decomp_bar,
    snp_coverage_bar,
    trait_limitations_badges,
    safe_float,
    risk_color,
    risk_badge,
    build_radar_chart_js,
    trust_tier_legend,
    portability_banner,
    evidence_letter,
    evidence_badge,
    trait_anchor_id,
    build_top_findings,
    build_summary_cards,
    build_clinvar_section,
    UI,
)


class TestSafeFloat:
    def test_valid_float(self):
        assert safe_float("3.14") == 3.14

    def test_integer(self):
        assert safe_float(42) == 42.0

    def test_none(self):
        assert safe_float(None) == 0.0

    def test_invalid_string(self):
        assert safe_float("abc") == 0.0

    def test_default_value(self):
        assert safe_float("abc", -1.0) == -1.0


class TestRiskColor:
    def test_high_positive(self):
        assert risk_color(2.5) == "#e74c3c"

    def test_high_negative(self):
        assert risk_color(-2.1) == "#e74c3c"

    def test_medium(self):
        assert risk_color(1.5) == "#f39c12"

    def test_low(self):
        assert risk_color(0.3) == "#27ae60"

    def test_inverted_never_returns_red(self):
        """Polarity-inverted traits (e.g. Morning chronotype, Cognitive
        function - risk_category=='high' means favorable) never represent
        danger even at their unfavorable end, so red must never appear."""
        for z in (-5.0, -2.0, -0.1, 0.0, 0.1, 2.0, 5.0):
            assert risk_color(z, inverted=True) != "#e74c3c"

    def test_inverted_large_positive_z_is_green(self):
        """A strongly favorable result (e.g. high z on Morning chronotype)
        must render green, not the red a large |z| would get by default."""
        assert risk_color(3.0, inverted=True) == "#27ae60"
        assert risk_color(3.0, inverted=False) == "#e74c3c"  # sanity: default path is unaffected

    def test_inverted_negative_z_is_amber_not_red(self):
        assert risk_color(-3.0, inverted=True) == "#f39c12"


class TestRiskBadge:
    def test_default_high_says_risk(self):
        html = risk_badge("high", UI["en"])
        assert UI["en"]["risk_high"] in html
        assert "#c0392b" in html  # red text

    def test_inverted_high_says_favorable_not_risk(self):
        html = risk_badge("high", UI["en"], inverted=True)
        assert "RISK" not in html.upper()
        assert UI["en"]["favorable_high"] in html
        assert "#1e8449" in html  # green text, same as the non-inverted "low" (good) styling

    def test_inverted_low_is_not_red(self):
        """The unfavorable end of an inverted trait still isn't red - these
        traits don't represent danger, just 'less of a good thing'."""
        html = risk_badge("low", UI["en"], inverted=True)
        assert "#c0392b" not in html
        assert UI["en"]["favorable_low"] in html

    def test_inverted_medium_is_neutral(self):
        html = risk_badge("medium", UI["en"], inverted=True)
        assert UI["en"]["favorable_medium"] in html

    def test_bilingual_es_inverted(self):
        html = risk_badge("high", UI["es"], inverted=True)
        assert UI["es"]["favorable_high"] in html
        assert "RIESGO" not in html.upper()


class TestComputePerTraitConfidence:
    """Test the composite confidence score computation."""

    def test_perfect_inputs(self):
        """100% SNPs, R²=0.99, slope=1.0, evidence A, low uncertainty → high confidence."""
        trait_entry = {
            "n_snps_used": 4,
            "n_snps_total": 4,
            "uncertainty_score": 0.2,
        }
        cal_entry = {
            "r_squared": 0.99,
            "calibration_slope": 1.0,
            "is_well_calibrated": True,
        }
        evidence_scores = [100, 100, 100, 100]  # All A-level
        score = compute_per_trait_confidence(trait_entry, cal_entry, None, evidence_scores)
        # snp: 100, cal: 100, ev: 100, uncert: 100 = avg 100
        assert score >= 90

    def test_inverted_calibration(self):
        """Inverted calibration slope should drastically lower confidence."""
        trait_entry = {
            "n_snps_used": 2,
            "n_snps_total": 4,
            "uncertainty_score": 1.0,
        }
        cal_entry = {
            "r_squared": 0.0,
            "calibration_slope": -1.02,
            "is_well_calibrated": False,
        }
        evidence_scores = [50, 50]
        score = compute_per_trait_confidence(trait_entry, cal_entry, None, evidence_scores)
        # snp: 50, cal: 0 (inverted), ev: 50, uncert: 0 (saturated) = avg 25
        assert score < 40

    def test_low_snp_coverage(self):
        """Low SNP coverage reduces confidence."""
        trait_entry = {
            "n_snps_used": 1,
            "n_snps_total": 4,
            "uncertainty_score": 0.6,
        }
        cal_entry = {
            "r_squared": 0.85,
            "calibration_slope": 1.05,
            "is_well_calibrated": True,
        }
        evidence_scores = [75]
        score = compute_per_trait_confidence(trait_entry, cal_entry, None, evidence_scores)
        # snp: 25, cal: 100, ev: 75, uncert: 80 = avg 70
        assert 50 <= score <= 85

    def test_no_calibration_data(self):
        """Missing calibration data defaults to 50 for cal_score."""
        trait_entry = {
            "n_snps_used": 2,
            "n_snps_total": 2,
            "uncertainty_score": 0.5,
        }
        score = compute_per_trait_confidence(trait_entry, None, None, None)
        # snp: 100, cal: 50, ev: 50, uncert: 100 = avg 75
        assert 60 <= score <= 90

    def test_saturated_uncertainty(self):
        """Saturated uncertainty (1.0) gives 0 for that component."""
        trait_entry = {
            "n_snps_used": 4,
            "n_snps_total": 4,
            "uncertainty_score": 1.0,
        }
        cal_entry = {
            "r_squared": 0.9,
            "calibration_slope": 0.95,
            "is_well_calibrated": True,
        }
        evidence_scores = [100, 100, 100, 100]
        score = compute_per_trait_confidence(trait_entry, cal_entry, None, evidence_scores)
        # snp: 100, cal: 100, ev: 100, uncert: 0 = avg 75
        assert 70 <= score <= 80

    def test_score_bounded_0_to_100(self):
        """Confidence score should always be in [0, 100]."""
        for snps in [(0, 4), (1, 4), (2, 4), (4, 4)]:
            for slope in [-1.0, 0.5, 1.0, 1.5]:
                for r2 in [0.0, 0.5, 0.9]:
                    for uncertainty in [0.0, 0.5, 1.0]:
                        trait_entry = {
                            "n_snps_used": snps[0],
                            "n_snps_total": snps[1],
                            "uncertainty_score": uncertainty,
                        }
                        cal_entry = {
                            "r_squared": r2,
                            "calibration_slope": slope,
                        }
                        score = compute_per_trait_confidence(trait_entry, cal_entry, None, [50])
                        assert 0 <= score <= 100, (
                            f"Score {score} out of bounds for snps={snps}, "
                            f"slope={slope}, r2={r2}, uncertainty={uncertainty}"
                        )


class TestConfidenceStars:
    """Test the star rating display."""

    def test_high_score(self):
        html = confidence_stars(100)
        assert "100%" in html
        assert "★" in html
        assert "#27ae60" in html  # Green

    def test_medium_score(self):
        html = confidence_stars(60)
        assert "60%" in html
        assert "#f39c12" in html  # Amber

    def test_low_score(self):
        html = confidence_stars(25)
        assert "25%" in html
        assert "#e74c3c" in html  # Red

    def test_zero_score(self):
        html = confidence_stars(0)
        assert "0%" in html
        assert "☆" in html  # Empty stars

    def test_all_stars_vary(self):
        """Star count should vary with score."""
        high = confidence_stars(95)
        low = confidence_stars(10)
        assert high.count("★") > low.count("★")


class TestCalibrationFlag:
    """Test the calibration quality badge."""

    def test_good_calibration(self):
        cal = {"r_squared": 0.95, "calibration_slope": 1.02}
        html = calibration_flag(cal)
        assert "GOOD" in html

    def test_fair_calibration(self):
        cal = {"r_squared": 0.6, "calibration_slope": 1.2}
        html = calibration_flag(cal)
        assert "FAIR" in html

    def test_inverted_calibration(self):
        cal = {"r_squared": 0.0, "calibration_slope": -1.02}
        html = calibration_flag(cal)
        assert "INVERTED" in html

    def test_poor_calibration(self):
        cal = {"r_squared": 0.3, "calibration_slope": 0.4}
        html = calibration_flag(cal)
        assert "POOR" in html

    def test_no_calibration_data(self):
        html = calibration_flag(None)
        assert "N/A" in html

    def test_slope_near_boundary(self):
        """Slope at 0.85 should be GOOD if R² is high."""
        cal = {"r_squared": 0.9, "calibration_slope": 0.85}
        html = calibration_flag(cal)
        assert "GOOD" in html


class TestTrustTier:
    """Test the trust tier classification."""

    def test_tier1_high_trust(self):
        """Well-calibrated trait with good confidence → TIER 1."""
        cal = {"r_squared": 0.9, "calibration_slope": 1.0, "is_well_calibrated": True}
        tier = trust_tier(confidence_score=72, cal_entry=cal, snp_ratio=0.5, uncertainty=0.9)
        assert tier == "TIER 1"

    def test_tier2_moderate(self):
        """Acceptable calibration, moderate confidence → TIER 2."""
        cal = {"r_squared": 0.6, "calibration_slope": 1.1, "is_well_calibrated": False}
        tier = trust_tier(confidence_score=55, cal_entry=cal, snp_ratio=0.5, uncertainty=0.9)
        assert tier == "TIER 2"

    def test_tier3_inverted(self):
        """Inverted calibration → TIER 3."""
        cal = {"r_squared": 0.0, "calibration_slope": -1.0, "is_well_calibrated": False}
        tier = trust_tier(confidence_score=30, cal_entry=cal, snp_ratio=0.5, uncertainty=1.0)
        assert tier == "TIER 3"

    def test_tier3_low_confidence(self):
        """Very low confidence → TIER 3 even with good calibration."""
        cal = {"r_squared": 0.9, "calibration_slope": 1.0, "is_well_calibrated": True}
        tier = trust_tier(confidence_score=30, cal_entry=cal, snp_ratio=0.5, uncertainty=1.0)
        assert tier == "TIER 3"

    def test_tier3_low_snp_coverage(self):
        """Low SNP coverage → TIER 3."""
        cal = {"r_squared": 0.9, "calibration_slope": 1.0, "is_well_calibrated": True}
        tier = trust_tier(confidence_score=70, cal_entry=cal, snp_ratio=0.3, uncertainty=0.5)
        assert tier == "TIER 3"


class TestTrustBadge:
    """Test the trust tier badge HTML."""

    def test_tier1_badge(self):
        html = trust_badge("TIER 1")
        assert "TIER 1" in html
        assert "#27ae60" in html  # Green

    def test_tier2_badge(self):
        html = trust_badge("TIER 2")
        assert "TIER 2" in html
        assert "#f39c12" in html  # Amber

    def test_tier3_badge(self):
        html = trust_badge("TIER 3")
        assert "TIER 3" in html
        assert "#e74c3c" in html  # Red


class TestMiniDecompBar:
    """Test the 3-layer uncertainty decomposition bar."""

    def test_effect_dominates(self):
        decomp = {"genotype_fraction": 0.1, "ancestry_fraction": 0.1, "effect_fraction": 0.8}
        html = mini_decomp_bar(decomp)
        assert "Gen" in html
        assert "Eff" in html

    def test_genotype_dominates(self):
        decomp = {"genotype_fraction": 0.7, "ancestry_fraction": 0.1, "effect_fraction": 0.2}
        html = mini_decomp_bar(decomp)
        assert "Gen" in html

    def test_no_decomp_data(self):
        html = mini_decomp_bar(None)
        assert "—" in html

    def test_empty_decomp(self):
        html = mini_decomp_bar({})
        assert "—" in html  # Empty dict is falsy, returns dash


class TestSnpCoverageBar:
    """Test the SNP coverage visual bar."""

    def test_full_coverage(self):
        html = snp_coverage_bar(4, 4)
        assert "4/4" in html
        assert "#27ae60" in html  # Green for 100%

    def test_half_coverage(self):
        html = snp_coverage_bar(1, 2)
        assert "1/2" in html
        assert "#f39c12" in html  # Amber for 50%

    def test_low_coverage(self):
        html = snp_coverage_bar(1, 4)
        assert "1/4" in html
        assert "#e74c3c" in html  # Red for 25%

    def test_zero_total(self):
        html = snp_coverage_bar(0, 0)
        assert "0/0" in html


class TestTraitLimitationsBadges:
    """Test the per-trait limitation badges."""

    def test_no_issues(self):
        trait_entry = {"n_snps_used": 4, "n_snps_total": 4, "uncertainty_score": 0.3}
        cal_entry = {"r_squared": 0.9, "calibration_slope": 1.0}
        html = trait_limitations_badges(trait_entry, cal_entry)
        assert "—" in html  # No issues

    def test_low_snps(self):
        trait_entry = {"n_snps_used": 1, "n_snps_total": 4, "uncertainty_score": 0.3}
        cal_entry = {"r_squared": 0.9, "calibration_slope": 1.0}
        html = trait_limitations_badges(trait_entry, cal_entry)
        assert "Low SNPs" in html

    def test_inverted_calibration_issue(self):
        trait_entry = {"n_snps_used": 3, "n_snps_total": 4, "uncertainty_score": 0.3}
        cal_entry = {"r_squared": 0.0, "calibration_slope": -1.0}
        html = trait_limitations_badges(trait_entry, cal_entry)
        assert "Inverted" in html

    def test_high_uncertainty(self):
        trait_entry = {"n_snps_used": 3, "n_snps_total": 4, "uncertainty_score": 0.9}
        cal_entry = {"r_squared": 0.9, "calibration_slope": 1.0}
        html = trait_limitations_badges(trait_entry, cal_entry)
        assert "uncert" in html.lower()

    def test_multiple_issues(self):
        trait_entry = {"n_snps_used": 1, "n_snps_total": 4, "uncertainty_score": 0.9}
        cal_entry = {"r_squared": 0.0, "calibration_slope": -1.0}
        html = trait_limitations_badges(trait_entry, cal_entry)
        assert "Low SNPs" in html
        assert "Inverted" in html
        assert "uncert" in html.lower()


# ═══════════════════════════════════════════════════════════════════════════════
# RADAR CHART TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestRadarChart:
    """Tests for the Chart.js interactive radar chart."""

    def _make_entries(self, n=9):
        traits = ["Bitter taste", "Blood pressure", "Caffeine", "Folate", "Glucose",
                  "Hair color", "Lactose", "Lipid", "Vitamin D"]
        return [{"trait": traits[i], "population_zscore": (i - 4) * 0.5,
                 "raw_score": (i - 4) * 0.5, "risk_category": "high" if abs(i - 4) > 2 else "medium",
                 "n_snps_used": 2, "n_snps_total": 4, "uncertainty_score": 0.5,
                 "population_percentile": 50.0}
                for i in range(n)]

    def test_renders_canvas(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "<canvas" in html
        assert "radarChart" in html

    def test_has_chart_js_config(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "new Chart(" in html
        assert "type: 'radar'" in html
        assert "<script>" in html

    def test_has_data_embedded(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "zScores" in html
        assert "fullNames" in html
        assert "tierLabels" in html

    def test_fallback_too_few_traits(self):
        entries = self._make_entries(2)
        html = build_radar_chart_js(entries, {})
        assert "Insufficient data" in html
        assert "<canvas" not in html

    def test_onclick_navigation(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "onClick" in html
        assert "scrollIntoView" in html

    def test_tooltips(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "Percentile" in html
        assert "Trust" in html

    def test_animation(self):
        entries = self._make_entries(9)
        html = build_radar_chart_js(entries, {})
        assert "animation" in html
        assert "duration" in html


# ═══════════════════════════════════════════════════════════════════════════════
# TRUST TIER LEGEND TEST
# ═══════════════════════════════════════════════════════════════════════════════

class TestTrustTierLegend:
    def test_legend_contains_all_tiers(self):
        html = trust_tier_legend()
        assert "T1" in html
        assert "T2" in html
        assert "T3" in html
        assert "High Trust" in html
        assert "Moderate" in html
        assert "Low Trust" in html


# ═══════════════════════════════════════════════════════════════════════════════
# PORTABILITY BANNER TEST
# ═══════════════════════════════════════════════════════════════════════════════

class TestPortabilityBanner:
    def test_shows_afr_bias(self):
        data = {"global_bias_index": 0.182, "most_biased": "AFR",
                "populations": [{"population": "AFR", "status": "LIMITED_PORTABILITY"}]}
        html = portability_banner(data)
        assert "AFR" in html

    def test_no_data_empty(self):
        assert portability_banner({}) == ""
        assert portability_banner(None) == ""


# ═══════════════════════════════════════════════════════════════════════════════
# TOP FINDINGS TEST (IMPROVEMENT_PLAN.md 1.1)
# ═══════════════════════════════════════════════════════════════════════════════

class TestEvidenceLetter:
    def test_a_grade(self):
        assert evidence_letter(100) == "A"

    def test_d_grade(self):
        assert evidence_letter(10) == "D"

    def test_badge_contains_letter(self):
        assert "Evidence B" in evidence_badge(75)


class TestTraitAnchorId:
    def test_slugifies_spaces_and_case(self):
        assert trait_anchor_id("Lactose Intolerance") == "trait-lactose-intolerance"

    def test_strips_non_alphanumeric(self):
        assert trait_anchor_id("Vitamin D (binding)") == "trait-vitamin-d-binding"


def _entry(trait, z, pctl, risk, n_used=2, n_total=3):
    return {"trait": trait, "population_zscore": z, "population_percentile": pctl,
            "risk_category": risk, "n_snps_used": n_used, "n_snps_total": n_total}


class TestBuildTopFindings:
    def test_ranks_by_priority_and_excludes_unscored_traits(self):
        entries = [
            _entry("High evidence strong signal", z=3.0, pctl=99, risk="high"),
            _entry("Weak signal", z=0.2, pctl=55, risk="medium"),
            _entry("Zero coverage trait", z=5.0, pctl=99, risk="high", n_used=0, n_total=3),
        ]
        evidence_lookup = {"high evidence strong signal": [100, 100],  # A
                            "weak signal": [25, 25]}                   # D
        html = build_top_findings(entries, UI["en"], evidence_lookup=evidence_lookup)

        assert "Zero coverage trait" not in html  # n_used=0 must be excluded
        first_idx = html.find("High evidence strong signal")
        second_idx = html.find("Weak signal")
        assert first_idx != -1 and second_idx != -1
        assert first_idx < second_idx  # higher |z| x evidence x confidence ranks first

    def test_anchor_link_targets_the_prs_table_row(self):
        entries = [_entry("Lactose intolerance", z=2.0, pctl=95, risk="high")]
        html = build_top_findings(entries, UI["en"],
                                   evidence_lookup={"lactose intolerance": [100]})
        assert f'href="#{trait_anchor_id("Lactose intolerance")}"' in html

    def test_no_fabricated_recommendation_when_none_curated(self):
        entries = [_entry("Some trait", z=2.0, pctl=95, risk="high")]
        html = build_top_findings(entries, UI["en"])
        assert UI["en"]["top_findings_action_fallback"] in html

    def test_uses_curated_recommendation_when_present_and_high_risk(self):
        entry = _entry("Some trait", z=2.0, pctl=95, risk="high")
        reco_lookup = {"some trait": {"recommendation_en": "Specific curated action text",
                                       "recommendation_es": "Texto de acción curado específico"}}
        html = build_top_findings([entry], UI["en"], recommendation_lookup=reco_lookup)
        assert "Specific curated action text" in html
        assert UI["en"]["top_findings_action_fallback"] not in html

    def test_curated_recommendation_ignored_when_risk_not_high(self):
        """A curated recommendation is written for the elevated-risk
        direction - it must not be shown for a medium/low (average or
        protective) finding just because the trait happens to be curated."""
        entry = _entry("Some trait", z=-2.0, pctl=5, risk="low")
        reco_lookup = {"some trait": {"recommendation_en": "Specific curated action text"}}
        html = build_top_findings([entry], UI["en"], recommendation_lookup=reco_lookup)
        assert "Specific curated action text" not in html
        assert UI["en"]["top_findings_action_fallback"] in html

    def test_empty_when_nothing_scored(self):
        html = build_top_findings([], UI["en"])
        assert UI["en"]["top_findings_empty"] in html

    def test_bilingual_es(self):
        entries = [_entry("Lactose intolerance", z=2.0, pctl=95, risk="high")]
        html = build_top_findings(entries, UI["es"],
                                   evidence_lookup={"lactose intolerance": [100]})
        assert UI["es"]["top_findings_action_fallback"] in html

    def test_polarity_inverted_trait_gets_favorable_badge_not_risk_badge(self):
        """The bug this project shipped for months: a favorable finding
        (e.g. Morning chronotype, risk_category=='high' meaning genetically
        morning-type/protective) rendered as a red 'HIGHER RISK' badge.
        Passing the trait in polarity_inverted must suppress that."""
        entries = [_entry("Morning chronotype (early bird)", z=2.5, pctl=98, risk="high")]
        html = build_top_findings(
            entries, UI["en"],
            evidence_lookup={"morning chronotype (early bird)": [100]},
            polarity_inverted={"morning chronotype (early bird)"},
        )
        assert UI["en"]["risk_high"] not in html
        assert UI["en"]["favorable_high"] in html

    def test_trait_not_in_polarity_inverted_set_is_unaffected(self):
        entries = [_entry("Lactose intolerance", z=2.5, pctl=98, risk="high")]
        html = build_top_findings(
            entries, UI["en"], evidence_lookup={"lactose intolerance": [100]},
            polarity_inverted={"morning chronotype (early bird)"},
        )
        assert UI["en"]["risk_high"] in html


class TestBuildSummaryCards:
    @staticmethod
    def _prs_result(entries):
        return {"prs_entries": entries}

    def test_high_risk_trait_counts_as_high(self):
        html = build_summary_cards(
            self._prs_result([_entry("Lactose intolerance", z=2.0, pctl=95, risk="high")]),
            ancestry={}, integrity={}, validation={}, ui=UI["en"],
        )
        assert ">1<" in html.split(UI["en"]["risk_high"])[0][-60:]

    def test_polarity_inverted_high_counts_as_favorable_not_risky(self):
        """A trait where 'high' is a good result must not inflate the
        executive summary's 'HIGHER RISK' count - that count is meant to
        flag traits worth a closer look for elevated risk, and a favorable
        finding is the opposite of that."""
        entries = [_entry("Morning chronotype (early bird)", z=2.0, pctl=95, risk="high")]
        html_default = build_summary_cards(
            self._prs_result(entries), ancestry={}, integrity={}, validation={}, ui=UI["en"],
        )
        html_inverted = build_summary_cards(
            self._prs_result(entries), ancestry={}, integrity={}, validation={}, ui=UI["en"],
            polarity_inverted={"morning chronotype (early bird)"},
        )
        # Without the polarity fix, this "high" trait is counted under risk_high.
        assert ">1<" in html_default.split(UI["en"]["risk_high"])[0][-60:]
        # With it, the same entry counts under risk_low instead.
        assert ">1<" in html_inverted.split(UI["en"]["risk_low"])[0][-60:]

    def test_medium_is_unaffected_by_polarity(self):
        entries = [_entry("Morning chronotype (early bird)", z=0.1, pctl=52, risk="medium")]
        html = build_summary_cards(
            self._prs_result(entries), ancestry={}, integrity={}, validation={}, ui=UI["en"],
            polarity_inverted={"morning chronotype (early bird)"},
        )
        assert ">1<" in html.split(UI["en"]["risk_medium"])[0][-60:]


# ═══════════════════════════════════════════════════════════════════════════════
# CURATED RECOMMENDATIONS DATA FILE TEST (IMPROVEMENT_PLAN.md 1.4)
# ═══════════════════════════════════════════════════════════════════════════════

class TestTraitRecommendationsData:
    """Regression coverage for data/trait_recommendations.json: every curated
    trait key must still exist as a real trait_category in the SNP panel, and
    every entry must have both languages + a reference. Catches the panel
    being re-curated (trait renamed/removed) out from under this file, and
    catches a future addition skipping the citation this session insisted on."""

    @staticmethod
    def _load():
        import json
        path = os.path.join(os.path.dirname(__file__), "..",
            "prs_research_pipeline", "data", "trait_recommendations.json")
        with open(path) as fh:
            return json.load(fh)

    @staticmethod
    def _real_trait_categories():
        import pandas as pd
        csv_path = os.path.join(os.path.dirname(__file__), "..",
            "prs_research_pipeline", "data", "snp_database_annotated.csv")
        return set(pd.read_csv(csv_path)["trait_category"].dropna().unique())

    def test_every_curated_trait_exists_in_the_panel(self):
        data = self._load()
        real_traits = self._real_trait_categories()
        for trait in data:
            if trait.startswith("_"):
                continue
            assert trait in real_traits, f"{trait!r} is not a real trait_category in the panel"

    def test_every_entry_has_both_languages_evidence_and_reference(self):
        data = self._load()
        for trait, entry in data.items():
            if trait.startswith("_"):
                continue
            for field in ("recommendation_en", "recommendation_es", "evidence_level", "reference"):
                assert entry.get(field), f"{trait!r} is missing {field!r}"

    def test_polarity_inverted_traits_exist_in_the_panel(self):
        """_polarity_inverted (consumed by risk_badge/risk_color for traits where
        'high' means favorable, not risky - IMPROVEMENT_PLAN.md follow-up fixed
        2026-08-06) must only name real trait_categories, same invariant as the
        curated recommendations above."""
        data = self._load()
        real_traits = self._real_trait_categories()
        inverted = data.get("_polarity_inverted", [])
        assert inverted, "_polarity_inverted should not be empty - Morning chronotype/Cognitive function are known cases"
        for trait in inverted:
            assert trait in real_traits, f"{trait!r} in _polarity_inverted is not a real trait_category in the panel"


# ═══════════════════════════════════════════════════════════════════════════════
# CLINVAR / MEDGEN LINK-OUT (IMPROVEMENT_PLAN.md 1.3)
# ═══════════════════════════════════════════════════════════════════════════════

class TestClinvarMedgenLink:
    """medgen_enrich.py now persists medgen_cui onto each variant (previously
    computed then discarded); the report should turn that into a real link
    to the source MedGen record instead of just showing unlinked prose."""

    def _clinvar_data(self, medgen_cui="C0018995"):
        variant = {
            "rsid": "rs1800562", "chrom": "6", "pos": 26093141,
            "genes": ["HFE"], "disease_name": "Hereditary hemochromatosis",
            "clinical_significance": "Pathogenic", "confidence_tier": "high",
            "review_status": "reviewed by expert panel",
            "disease_description": "A disorder of iron metabolism.",
        }
        if medgen_cui:
            variant["medgen_cui"] = medgen_cui
        return {
            "pathogenic_variants": [variant],
            "pathogenic_variant_summary": {"high_confidence_count": 1, "by_confidence_tier": {}},
            "metadata": {},
        }

    def test_renders_medgen_link_when_cui_present(self):
        html = build_clinvar_section(self._clinvar_data(medgen_cui="C0018995"), UI["en"])
        assert "ncbi.nlm.nih.gov/medgen/C0018995" in html

    def test_no_broken_link_when_cui_absent(self):
        html = build_clinvar_section(self._clinvar_data(medgen_cui=""), UI["en"])
        assert "ncbi.nlm.nih.gov/medgen/" not in html


class TestClinvarBodySystemGrouping:
    """medgen_enrich.py persists body_system on each variant (IMPROVEMENT_PLAN.md
    1.3 - "agrupar por sistema/organo"); the report should render a system
    header per group instead of one flat table in raw file order."""

    def _clinvar_data(self, variants):
        return {
            "pathogenic_variants": variants,
            "pathogenic_variant_summary": {"high_confidence_count": len(variants), "by_confidence_tier": {}},
            "metadata": {},
        }

    def _variant(self, rsid, body_system, confidence_tier="high"):
        return {
            "rsid": rsid, "chrom": "1", "pos": 1, "genes": ["X"],
            "disease_name": "not_provided", "clinical_significance": "Pathogenic",
            "confidence_tier": confidence_tier, "body_system": body_system,
        }

    def test_renders_a_header_per_distinct_system(self):
        html = build_clinvar_section(self._clinvar_data([
            self._variant("rs1", "oncology"),
            self._variant("rs2", "cardiovascular"),
        ]), UI["en"])
        assert "Oncology" in html
        assert "Cardiovascular" in html

    def test_larger_group_renders_before_smaller_group(self):
        html = build_clinvar_section(self._clinvar_data([
            self._variant("rs1", "cardiovascular"),
            self._variant("rs2", "oncology"),
            self._variant("rs3", "oncology"),
        ]), UI["en"])
        assert html.index("Oncology") < html.index("Cardiovascular")

    def test_missing_body_system_falls_back_to_other_not_a_crash(self):
        variant = self._variant("rs1", body_system=None)
        del variant["body_system"]
        html = build_clinvar_section(self._clinvar_data([variant]), UI["en"])
        assert "Other" in html
