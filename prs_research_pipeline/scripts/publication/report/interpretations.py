"""
Pure interpretation helpers extracted from comprehensive_report.py
(IMPROVEMENT_PLAN.md 1.6, Phase 1). Every function here turns data into a
decision/badge/color/tier - moved verbatim, no logic changes.
"""

import json
import re
from pathlib import Path


def load_json(path: str, default=None):
    p = Path(path)
    if p.exists():
        with open(p) as fh:
            return json.load(fh)
    return default or {}

def safe_float(v, default=0.0):
    try: return float(v)
    except (ValueError, TypeError): return default


def trait_anchor_id(trait: str) -> str:
    """Slugify a trait name into an HTML id usable for in-page anchor links."""
    slug = re.sub(r"[^a-z0-9]+", "-", trait.lower()).strip("-")
    return f"trait-{slug}"

def risk_color(z_score, inverted=False):
    """Return CSS color based on z-score.

    `inverted=True` is for traits where a high score is a favorable
    finding, not elevated risk (e.g. Morning chronotype, Cognitive
    function - see _polarity_inverted in trait_recommendations.json).
    For those, color follows direction (green=favorable, amber=less
    favorable) rather than |z| magnitude, and never turns red - these
    traits don't represent danger even at their "low" end.
    """
    if inverted:
        return "#27ae60" if z_score >= 0 else "#f39c12"
    az = abs(z_score)
    if az >= 2.0: return "#e74c3c"
    if az >= 1.0: return "#f39c12"
    return "#27ae60"

def risk_badge(risk, ui, inverted=False):
    """Render the risk/direction badge.

    `inverted=True` swaps the badge to neutral "favorable/typical/less
    favorable" wording and colors instead of "risk" language, for
    traits where risk_category=="high" means a good result (see
    risk_color docstring).
    """
    if inverted:
        labels = {"high": ui["favorable_high"], "medium": ui["favorable_medium"], "low": ui["favorable_low"]}
        colors = {"high": ("#d5f5e3", "#1e8449"), "medium": ("#fdebd0", "#b7950b"), "low": ("#fdebd0", "#b7950b")}
    else:
        labels = {"high": ui["risk_high"], "medium": ui["risk_medium"], "low": ui["risk_low"]}
        colors = {"high": ("#fadbd8", "#c0392b"), "medium": ("#fdebd0", "#b7950b"), "low": ("#d5f5e3", "#1e8449")}
    bg, fg = colors.get(risk, colors["medium"])
    return f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.75rem;font-weight:700">{labels.get(risk, risk)}</span>'

def risk_bar(pct, z_score, inverted=False):
    """Render a horizontal risk bar."""
    color = risk_color(z_score, inverted=inverted)
    return (
        f'<div style="display:flex;align-items:center;gap:6px;margin:4px 0">'
        f'<span style="font-size:0.65rem;color:#27ae60">Low</span>'
        f'<div style="flex:1;height:8px;background:#e9ecef;border-radius:4px;overflow:hidden">'
        f'<div style="width:{pct}%;height:100%;background:{color};border-radius:4px"></div>'
        f'</div>'
        f'<span style="font-size:0.65rem;color:#e74c3c">High</span>'
        f'</div>'
    )


def compute_per_trait_confidence(trait_entry, cal_entry, uncert_entry, evidence_scores=None):
    """Compute a composite 0-100 confidence score for a trait.

    Components (equal weight):
      - snp_coverage_score: fraction of expected SNPs found
      - calibration_score: based on R² and slope correctness
      - evidence_score: average evidence level (A=100, B=75, C=50, D=25)
      - uncertainty_score: inverted (1.0 - uncertainty_score) scaled to 0-100
    """
    n_used = trait_entry.get("n_snps_used", 0)
    n_total = trait_entry.get("n_snps_total", 1)
    snp_ratio = n_used / max(n_total, 1)

    # SNP coverage: linear scaling
    snp_score = snp_ratio * 100

    # Calibration quality: based on R² and slope
    if cal_entry:
        r2 = safe_float(cal_entry.get("r_squared", 0))
        slope = safe_float(cal_entry.get("calibration_slope", 1.0))
        if slope < 0:
            cal_score = 0  # Inverted calibration = zero confidence
        elif r2 >= 0.8 and 0.85 <= slope <= 1.15:
            cal_score = 100
        elif r2 >= 0.5 and 0.5 <= slope <= 1.5:
            cal_score = 60
        else:
            cal_score = max(0, r2 * 100)
    else:
        cal_score = 50  # No calibration data available

    # Evidence level: average across all SNPs
    if evidence_scores:
        ev_score = sum(evidence_scores) / max(len(evidence_scores), 1)
    else:
        ev_score = 50  # Unknown evidence

    # Uncertainty: inverted
    uncertainty = safe_float(trait_entry.get("uncertainty_score", 1.0))
    if uncertainty >= 1.0:
        uncert_score = 0
    elif uncertainty <= 0.5:
        uncert_score = 100
    else:
        uncert_score = (1.0 - uncertainty) * 200

    # Composite with equal weights
    confidence = 0.25 * snp_score + 0.25 * cal_score + 0.25 * ev_score + 0.25 * uncert_score
    return round(max(0, min(100, confidence)))


def confidence_stars(score):
    """Render 0-5 star confidence indicator."""
    stars = score / 20  # 0-100 -> 0-5
    full = int(stars)
    half = 1 if (stars - full) >= 0.5 else 0
    empty = 5 - full - half
    color = "#27ae60" if score >= 75 else ("#f39c12" if score >= 50 else "#e74c3c")
    # Half star: a full star (★) clipped to 50% width over an empty star (☆)
    # background, rather than the U+2BE8 "left half black star" glyph - that
    # codepoint has essentially no font coverage and renders as a broken
    # tofu/fallback box in most browsers, despite looking fine in an editor.
    half_star = (
        '<span style="position:relative;display:inline-block;width:1em">☆'
        '<span style="position:absolute;left:0;top:0;width:0.5em;overflow:hidden;'
        'display:inline-block">★</span></span>'
    ) if half else ""
    return (
        f'<div style="display:flex;align-items:center;gap:4px;min-width:110px">'
        f'<span style="color:{color};font-weight:700;font-size:0.85rem;min-width:2.5em">{score:.0f}%</span>'
        f'<span style="color:#f1c40f;font-size:0.75rem">{"★" * full}{half_star}{"☆" * empty}</span>'
        f'</div>'
    )


def calibration_flag(cal_entry):
    """Render calibration quality badge for a trait."""
    if not cal_entry:
        return '<span style="color:#95a5a6;font-size:0.7rem">N/A</span>'
    slope = safe_float(cal_entry.get("calibration_slope", 1.0))
    r2 = safe_float(cal_entry.get("r_squared", 0.0))
    if slope < 0:
        return (
            f'<span style="background:#fadbd8;color:#c0392b;padding:2px 6px;border-radius:3px;'
            f'font-size:0.65rem;font-weight:700" '
            f'title="Direction reversed! Slope={slope:.2f}">INVERTED</span>'
        )
    if r2 >= 0.8 and 0.85 <= slope <= 1.15:
        return (
            f'<span style="background:#d5f5e3;color:#1e8449;padding:2px 6px;border-radius:3px;'
            f'font-size:0.65rem;font-weight:700" '
            f'title="R²={r2:.3f}, slope={slope:.2f}">GOOD</span>'
        )
    if r2 >= 0.5:
        return (
            f'<span style="background:#fdebd0;color:#b7950b;padding:2px 6px;border-radius:3px;'
            f'font-size:0.65rem;font-weight:700" '
            f'title="R²={r2:.3f}, slope={slope:.2f}">FAIR</span>'
        )
    return (
        f'<span style="background:#fadbd8;color:#c0392b;padding:2px 6px;border-radius:3px;'
        f'font-size:0.65rem;font-weight:700" '
        f'title="R²={r2:.3f}, slope={slope:.2f}">POOR</span>'
    )


def trust_tier(confidence_score, cal_entry, snp_ratio, uncertainty):
    """Classify a trait into TIER 1 (High), TIER 2 (Moderate), or TIER 3 (Low) trust.

    Thresholds are calibrated for the current data reality:
    - Most traits have 50% SNP coverage (1-2 of 2-4 SNPs available)
    - Uncertainty is often saturated at 1.0 (effect SE dominates)
    - TIER 1: good calibration + reasonable confidence despite data limitations
    - TIER 2: acceptable calibration (not inverted) + moderate confidence
    - TIER 3: inverted calibration or very low confidence
    """
    # TIER 1: >=50% SNPs, good calibration, confidence >= 60
    if (snp_ratio >= 0.5 and
            cal_entry and cal_entry.get("is_well_calibrated", False) and
            safe_float(cal_entry.get("calibration_slope", 1.0)) >= 0 and
            confidence_score >= 60):
        return "TIER 1"
    # TIER 2: >=50% SNPs, NOT inverted calibration, confidence >= 40
    if (snp_ratio >= 0.5 and
            cal_entry and
            safe_float(cal_entry.get("calibration_slope", 0)) >= 0 and
            confidence_score >= 40):
        return "TIER 2"
    return "TIER 3"


def trust_badge(tier):
    """Render a trust tier badge."""
    styles = {
        "TIER 1": ("#27ae60", "#d5f5e3", "High Trust"),
        "TIER 2": ("#f39c12", "#fdebd0", "Moderate Trust"),
        "TIER 3": ("#e74c3c", "#fadbd8", "Low Trust"),
    }
    fg, bg, label = styles.get(tier, ("#95a5a6", "#eaecee", "Unknown"))
    return (
        f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;'
        f'font-size:0.65rem;font-weight:700;white-space:nowrap" '
        f'title="{tier}: {label}">{tier}</span>'
    )


def mini_decomp_bar(decomp):
    """Render a compact 3-layer stacked uncertainty bar."""
    if not decomp:
        return '<span style="color:#95a5a6;font-size:0.7rem">—</span>'
    g = max(safe_float(decomp.get("genotype_fraction", 0)) * 100, 1)
    a = max(safe_float(decomp.get("ancestry_fraction", 0)) * 100, 1)
    e = max(safe_float(decomp.get("effect_fraction", 0)) * 100, 1)
    # Dominant component indicator
    if e > 50:
        dot_color, dot_title = "#e74c3c", "Effect SE dominates uncertainty"
    elif g > 50:
        dot_color, dot_title = "#3498db", "Genotype quality drives uncertainty"
    else:
        dot_color, dot_title = "#f39c12", "Ancestry ambiguity drives uncertainty"

    return (
        f'<div style="display:flex;flex-direction:column;gap:2px;min-width:80px">'
        f'<div style="display:flex;height:8px;border-radius:3px;overflow:hidden;background:#e9ecef">'
        f'<div style="width:{g:.0f}%;background:#3498db" title="Genotype: {g:.0f}%"></div>'
        f'<div style="width:{a:.0f}%;background:#e74c3c" title="Ancestry: {a:.0f}%"></div>'
        f'<div style="width:{e:.0f}%;background:#f39c12" title="Effect: {e:.0f}%"></div>'
        f'</div>'
        f'<div style="display:flex;justify-content:space-between;font-size:0.55rem;color:#7f8c8d">'
        f'<span>Gen</span><span style="color:{dot_color};font-weight:700" title="{dot_title}">●</span><span>Eff</span>'
        f'</div>'
        f'</div>'
    )


def snp_coverage_bar(n_used, n_total):
    """Render a compact SNP coverage bar with ratio."""
    ratio = n_used / max(n_total, 1) * 100
    color = "#27ae60" if ratio >= 75 else ("#f39c12" if ratio >= 50 else "#e74c3c")
    return (
        f'<div style="display:flex;align-items:center;gap:4px;min-width:70px">'
        f'<div style="flex:1;height:6px;background:#e9ecef;border-radius:3px;overflow:hidden">'
        f'<div style="width:{ratio:.0f}%;height:100%;background:{color};border-radius:3px"></div>'
        f'</div>'
        f'<span style="font-size:0.7rem;font-weight:600;color:{color};white-space:nowrap">{n_used}/{n_total}</span>'
        f'</div>'
    )


def portability_banner(portability_data):
    """Build a portability warning banner for the PRS section."""
    if not portability_data:
        return ""
    global_bias = safe_float(portability_data.get("global_bias_index", 0))
    pops = portability_data.get("populations", [])

    pop_chips = ""
    for p in pops:
        pop = p.get("population", "")
        status = p.get("status", "")
        colors = {
            "GOOD_PORTABILITY": ("#d5f5e3", "#1e8449"),
            "MODERATE_PORTABILITY": ("#fdebd0", "#b7950b"),
            "LIMITED_PORTABILITY": ("#fadbd8", "#c0392b"),
        }
        bg, fg = colors.get(status, ("#eaecee", "#7f8c8d"))
        label = status.replace("_PORTABILITY", "").title() if status else "?"
        pop_chips += (
            f'<span style="background:{bg};color:{fg};padding:1px 6px;border-radius:3px;'
            f'font-size:0.65rem;font-weight:700">{pop}: {label}</span> '
        )

    return (
        f'<div style="background:#fef9e7;border:2px solid #f39c12;border-radius:8px;'
        f'padding:0.8rem 1.2rem;margin-bottom:1rem">'
        f'<div style="display:flex;align-items:flex-start;gap:8px">'
        f'<span style="font-size:1.2rem">⚠️</span>'
        f'<div>'
        f'<strong style="color:#b7950b">Population Portability Notice:</strong>'
        f'<p style="font-size:0.8rem;margin:4px 0 0;color:#7d6608">'
        f'PRS scores are calibrated for EUR populations. Cross-population portability is '
        f'<strong>LIMITED</strong> for non-EUR ancestries '
        f'(Global Bias Index: {global_bias:.3f}). '
        f'AFR shows the highest bias due to differing LD structure and allele frequencies. '
        f'Results should be interpreted with caution for non-European samples.'
        f'</p>'
        f'<div style="margin-top:6px;display:flex;gap:6px;flex-wrap:wrap">{pop_chips}</div>'
        f'</div></div></div>'
    )


def reference_coverage_banner(prs_result, lang="en"):
    """Prominent banner when Stage C fell back to the chr22-only 1000G reference
    (IMPROVEMENT_PLAN.md 0.2). PCA-based ancestry and population calibration are
    computed against that reference, so a chr22-only run means those components
    reflect a single chromosome, not the genome — worth surfacing before any
    other section, not just in the debug log."""
    coverage = prs_result.get("metadata", {}).get("reference_coverage", "genome_wide")
    if coverage != "chr22_only":
        return ""
    if lang == "es":
        title = "Aviso: cobertura de referencia parcial (solo chr22)"
        body = (
            "Este informe se generó sin la referencia 1000G genoma-completo instalada, así que "
            "la ancestría (PCA) y la calibración poblacional se calcularon usando <strong>solo "
            "el cromosoma 22</strong> como referencia, no el genoma completo. Los rasgos cuyos "
            "SNPs caen fuera de chr22 pueden faltar o estar mal calibrados. "
            "Ejecuta <code>scripts/setup/download_1000G_full.py</code> y vuelve a correr el "
            "pipeline para la cobertura completa."
        )
    else:
        title = "Notice: Partial Reference Coverage (chr22 only)"
        body = (
            "This report was generated without the genome-wide 1000G reference installed, so "
            "ancestry (PCA) and population calibration were computed using <strong>chromosome 22 "
            "only</strong> as the reference, not the full genome. Traits whose SNPs fall outside "
            "chr22 may be missing or miscalibrated. "
            "Run <code>scripts/setup/download_1000G_full.py</code> and re-run the pipeline for "
            "full coverage."
        )
    return (
        f'<div class="portability-banner" style="padding:0.8rem 1.2rem;margin-bottom:1rem">'
        f'<div style="display:flex;align-items:flex-start;gap:8px">'
        f'<span style="font-size:1.2rem">⚠️</span>'
        f'<div>'
        f'<strong style="color:#b7950b">{title}</strong>'
        f'<p style="font-size:0.8rem;margin:4px 0 0;color:#7d6608">{body}</p>'
        f'</div></div></div>'
    )


def trait_limitations_badges(trait_entry, cal_entry):
    """Generate per-trait limitation badges for the PRS table."""
    issues = []
    n_used = trait_entry.get("n_snps_used", 0)
    n_total = trait_entry.get("n_snps_total", 0)
    uncertainty = safe_float(trait_entry.get("uncertainty_score", 1.0))

    if n_total > 0 and n_used / n_total < 0.5:
        issues.append(("Low SNPs", "#fadbd8", "#c0392b"))
    if cal_entry and safe_float(cal_entry.get("calibration_slope", 1.0)) < 0:
        issues.append(("Inverted cal.", "#fdebd0", "#b7950b"))
    if uncertainty >= 0.8:
        issues.append(("High uncert.", "#fdebd0", "#b7950b"))

    if not issues:
        return '<span style="color:#27ae60;font-size:0.65rem">—</span>'

    badges = ""
    for text, bg, fg in issues:
        badges += (
            f'<span style="background:{bg};color:{fg};padding:1px 5px;border-radius:3px;'
            f'font-size:0.6rem;font-weight:700;margin-right:2px;white-space:nowrap">{text}</span>'
        )
    return badges


def trust_tier_legend():
    """Render a legend box explaining the trust tiers."""
    return (
        '<div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;'
        'padding:0.6rem 1rem;margin-bottom:0.8rem;font-size:0.72rem">'
        '<strong style="margin-right:8px">Trust Tiers:</strong>'
        '<span style="background:#d5f5e3;color:#1e8449;padding:2px 6px;border-radius:3px;'
        'font-weight:700;margin-right:4px">T1 High Trust</span> '
        'Good calibration + ≥50% SNPs + confidence ≥60% | '
        '<span style="background:#fdebd0;color:#b7950b;padding:2px 6px;border-radius:3px;'
        'font-weight:700;margin-right:4px">T2 Moderate</span> '
        '≥50% SNPs + not inverted + confidence ≥40% | '
        '<span style="background:#fadbd8;color:#c0392b;padding:2px 6px;border-radius:3px;'
        'font-weight:700;margin-right:4px">T3 Low Trust</span> '
        'Inverted calibration or very low confidence'
        '</div>'
    )


def evidence_letter(evidence_avg_score):
    """Map an averaged 0-100 evidence score (A=100,B=75,C=50,D=25) back to a letter."""
    if evidence_avg_score >= 90: return "A"
    if evidence_avg_score >= 65: return "B"
    if evidence_avg_score >= 40: return "C"
    return "D"

def evidence_badge(evidence_avg_score):
    letter = evidence_letter(evidence_avg_score)
    colors = {"A": ("#d5f5e3", "#1e8449"), "B": ("#d6eaf8", "#2874a6"),
              "C": ("#fdebd0", "#b7950b"), "D": ("#fadbd8", "#c0392b")}
    bg, fg = colors[letter]
    return f'<span style="background:{bg};color:{fg};padding:2px 7px;border-radius:4px;font-size:0.7rem;font-weight:700">Evidence {letter}</span>'
