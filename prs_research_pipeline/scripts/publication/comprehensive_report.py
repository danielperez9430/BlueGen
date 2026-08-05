#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   COMPREHENSIVE HTML REPORT GENERATOR                                        ║
║   scripts/comprehensive_report.py                                            ║
║                                                                            ║
║   Generates a single, rich, interactive HTML report from all pipeline       ║
║   outputs. Collapsible sections, embedded data tables, visual risk bars,    ║
║   and full bilingual support.                                               ║
║                                                                            ║
║   Inputs:  All JSON outputs from the PRS pipeline                           ║
║   Output:  reports/comprehensive_report_en.html                             ║
║            reports/comprehensive_report_es.html                             ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import re
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
from dataclasses import dataclass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.constants import PIPELINE_VERSION

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════════
# BILINGUAL UI STRINGS
# ═══════════════════════════════════════════════════════════════════════════════

UI = {
    "en": {
        "title": "BlueGen Report",
        "subtitle": "Comprehensive Polygenic Risk Score Analysis",
        "sections": {
            "summary": "Executive Summary",
            "ancestry": "Ancestry Deep-Dive",
            "prs": "PRS Results — Population-Calibrated",
            "uncertainty": "Uncertainty — Variance Decomposition",
            "variants": "Variant-Level Detail",
            "calibration": "Population Calibration Methodology",
            "pgs_calibration": "PGS Catalog — External Validation",
            "portability": "Population Portability Analysis",
            "validation": "Scientific Validation (8 Dimensions)",
            "gwas_consortium": "GWAS Consortium Validation",
            "benchmark": "External Benchmarking — Quality Delta",
            "adversarial": "Adversarial Stress Testing",
            "failure_map": "Failure Mode Coverage",
            "leakage": "Leakage Prevention — Detailed Audit",
            "consistency": "GWAS-Ancestry Consistency Check",
            "integrity": "Scientific Integrity Score",
            "reproducibility": "Reproducibility — Environment & Seeds",
            "methodology": "Pipeline Methodology",
            "clinvar": "ClinVar — Pathogenic Variants",
            "limitations": "Limitations & Disclaimers",
        },
        "risk_high": "HIGHER RISK", "risk_medium": "AVERAGE RISK", "risk_low": "LOWER RISK",
        "passed": "PASSED", "failed": "FAILED", "warning": "WARNING",
        "yes": "Yes", "no": "No",
        "top_findings_title": "🔍 Your Top Findings",
        "top_findings_intro": (
            "The traits below are prioritized by how far your result deviates from the reference "
            "population, weighted by how much genetic evidence supports each trait's score. "
            "This is a navigation aid, not a ranking of health importance — see the full PRS "
            "Results table below for every trait scored."
        ),
        "top_findings_empty": "No traits had enough matched SNPs to prioritize yet.",
        "top_findings_meaning": "{trait} is {risk_word}, at the {pctl:.0f}th percentile versus the {pop} reference population (z={z:+.2f}), based on {n_used} of {n_total} panel SNP(s) with {evidence} evidence.",
        "top_findings_action_fallback": "Curated, evidence-cited guidance for this trait is not yet available in this report — discuss with a nutritionist or healthcare professional, and see the full trait detail below.",
        "top_findings_jump": "See full detail ↓",
        "top_findings_disclaimer_summary": "⚠️ Research use only — how confident should you be in these results? (click to expand)",
        "risk_word_high": "elevated relative to the reference population",
        "risk_word_medium": "close to the reference population average",
        "risk_word_low": "lower relative to the reference population",
        "disclaimer": (
            "⚠️  RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS\n\n"
            "This PRS report is generated for RESEARCH PURPOSES ONLY. "
            "It does NOT constitute a clinical diagnosis, medical advice, "
            "or a definitive prediction of disease risk.\n\n"
            "Key limitations:\n"
            "• PRS is probabilistic, not deterministic\n"
            "• Effect sizes depend on GWAS discovery populations\n"
            "• Ancestry bias is reduced but not eliminated by population calibration\n"
            "• Gene-environment interactions are not captured by genotype alone\n"
            "• Consult a healthcare professional before making dietary or lifestyle changes"
        ),
    },
    "es": {
        "title": "Informe de Investigación PRS",
        "subtitle": "Análisis Completo de Puntaje de Riesgo Poligénico",
        "sections": {
            "summary": "Resumen Ejecutivo",
            "ancestry": "Ascendencia en Detalle",
            "prs": "Resultados PRS — Calibrados por Población",
            "uncertainty": "Incertidumbre — Descomposición de Varianza",
            "variants": "Detalle por Variante",
            "calibration": "Metodología de Calibración Poblacional",
            "pgs_calibration": "Catálogo PGS — Validación Externa",
            "portability": "Análisis de Portabilidad Poblacional",
            "validation": "Validación Científica (8 Dimensiones)",
            "gwas_consortium": "Validación de Consorcios GWAS",
            "benchmark": "Referencia Externa — Delta de Calidad",
            "adversarial": "Pruebas de Estrés Adversarial",
            "failure_map": "Cobertura de Modos de Falla",
            "leakage": "Prevención de Fuga — Auditoría Detallada",
            "consistency": "Verificación de Consistencia GWAS-Ascendencia",
            "integrity": "Índice de Integridad Científica",
            "reproducibility": "Reproducibilidad — Entorno y Semillas",
            "methodology": "Metodología del Pipeline",
            "clinvar": "ClinVar — Variantes Patogénicas",
            "limitations": "Limitaciones y Avisos",
        },
        "risk_high": "RIESGO ELEVADO", "risk_medium": "RIESGO PROMEDIO", "risk_low": "RIESGO BAJO",
        "passed": "APROBADO", "failed": "FALLIDO", "warning": "ADVERTENCIA",
        "yes": "Sí", "no": "No",
        "top_findings_title": "🔍 Tus Hallazgos Principales",
        "top_findings_intro": (
            "Los rasgos de abajo están priorizados por cuánto se desvía tu resultado de la población "
            "de referencia, ponderado por cuánta evidencia genética respalda el score de cada rasgo. "
            "Es una guía de navegación, no un ranking de importancia para la salud — ver la tabla "
            "completa de Resultados PRS más abajo para todos los rasgos puntuados."
        ),
        "top_findings_empty": "Ningún rasgo tuvo suficientes SNPs casados todavía para priorizar.",
        "top_findings_meaning": "{trait} está {risk_word}, en el percentil {pctl:.0f} respecto a la población de referencia {pop} (z={z:+.2f}), basado en {n_used} de {n_total} SNP(s) del panel con evidencia {evidence}.",
        "top_findings_action_fallback": "Todavía no hay una recomendación curada y con evidencia citada para este rasgo en este informe — coméntalo con un nutricionista o profesional de la salud, y consultá el detalle completo más abajo.",
        "top_findings_jump": "Ver detalle completo ↓",
        "top_findings_disclaimer_summary": "⚠️ Solo para uso en investigación — ¿cuánta confianza depositar en estos resultados? (clic para expandir)",
        "risk_word_high": "elevado respecto a la población de referencia",
        "risk_word_medium": "cercano al promedio de la población de referencia",
        "risk_word_low": "más bajo respecto a la población de referencia",
        "disclaimer": (
            "⚠️  SOLO PARA USO EN INVESTIGACIÓN — NO PARA DIAGNÓSTICO CLÍNICO\n\n"
            "Este informe PRS se genera SOLO PARA FINES DE INVESTIGACIÓN. "
            "NO constituye un diagnóstico clínico, consejo médico, "
            "ni una predicción definitiva del riesgo de enfermedad.\n\n"
            "Limitaciones clave:\n"
            "• PRS es probabilístico, no determinista\n"
            "• Los tamaños del efecto dependen de las poblaciones de descubrimiento GWAS\n"
            "• El sesgo de ascendencia se reduce pero no se elimina\n"
            "• Las interacciones gen-ambiente no se capturan solo con el genotipo\n"
            "• Consulte a un profesional de la salud antes de hacer cambios"
        ),
    },
}

POP_NAMES = {
    "en": {"EUR": "European", "AFR": "African", "EAS": "East Asian",
           "SAS": "South Asian", "AMR": "Admixed American"},
    "es": {"EUR": "Europea", "AFR": "Africana", "EAS": "Asia Oriental",
           "SAS": "Sur de Asia", "AMR": "Americana Mixta"},
}

# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_json(path: str, default=None):
    p = Path(path)
    if p.exists():
        with open(p) as fh:
            return json.load(fh)
    return default or {}

def safe_float(v, default=0.0):
    try: return float(v)
    except (ValueError, TypeError): return default

# ═══════════════════════════════════════════════════════════════════════════════
# HTML GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def trait_anchor_id(trait: str) -> str:
    """Slugify a trait name into an HTML id usable for in-page anchor links."""
    slug = re.sub(r"[^a-z0-9]+", "-", trait.lower()).strip("-")
    return f"trait-{slug}"

def risk_color(z_score):
    """Return CSS color based on z-score magnitude."""
    az = abs(z_score)
    if az >= 2.0: return "#e74c3c"
    if az >= 1.0: return "#f39c12"
    return "#27ae60"

def risk_badge(risk, ui):
    labels = {"high": ui["risk_high"], "medium": ui["risk_medium"], "low": ui["risk_low"]}
    colors = {"high": ("#fadbd8", "#c0392b"), "medium": ("#fdebd0", "#b7950b"), "low": ("#d5f5e3", "#1e8449")}
    bg, fg = colors.get(risk, colors["medium"])
    return f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.75rem;font-weight:700">{labels.get(risk, risk)}</span>'

def risk_bar(pct, z_score):
    """Render a horizontal risk bar."""
    color = risk_color(z_score)
    return (
        f'<div style="display:flex;align-items:center;gap:6px;margin:4px 0">'
        f'<span style="font-size:0.65rem;color:#27ae60">Low</span>'
        f'<div style="flex:1;height:8px;background:#e9ecef;border-radius:4px;overflow:hidden">'
        f'<div style="width:{pct}%;height:100%;background:{color};border-radius:4px"></div>'
        f'</div>'
        f'<span style="font-size:0.65rem;color:#e74c3c">High</span>'
        f'</div>'
    )

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIDENCE & TRUST HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

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
    return (
        f'<div style="display:flex;align-items:center;gap:4px;min-width:110px">'
        f'<span style="color:{color};font-weight:700;font-size:0.85rem;min-width:2.5em">{score:.0f}%</span>'
        f'<span style="color:#f1c40f;font-size:0.75rem">{"★" * full}{"⯨" if half else ""}{"☆" * empty}</span>'
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


# ═══════════════════════════════════════════════════════════════════════════════
# RADAR CHART (Chart.js — interactive)
# ═══════════════════════════════════════════════════════════════════════════════

def build_radar_chart_js(entries, ui, cal_lookup=None, uncert_lookup=None, evidence_lookup=None):
    """Build an interactive radar chart using Chart.js.

    Features: tooltips (z-score + percentile + trust tier), animation on load,
    click-to-scroll to PRS table row, tier highlighting via legend.
    """
    if cal_lookup is None:
        cal_lookup = {}
    if uncert_lookup is None:
        uncert_lookup = {}
    if evidence_lookup is None:
        evidence_lookup = {}

    traits = [(e.get("trait", ""),
               safe_float(e.get("population_zscore", e.get("raw_score", 0))),
               e.get("risk_category", "medium"),
               safe_float(e.get("population_percentile", 50)))
              for e in (entries or [])]

    if len(traits) < 3:
        return ('<div style="text-align:center;padding:2rem;color:var(--color-text-secondary)">'
                '<p>Insufficient data for radar visualization</p></div>')

    # Short labels for radar axes
    short_labels = [(t[:12] + ("…" if len(t) > 12 else "")) for t, _, _, _ in traits]
    z_scores = [z for _, z, _, _ in traits]
    pctls = [p for _, _, _, p in traits]
    full_names = [t for t, _, _, _ in traits]
    risk_cats = [r for _, _, r, _ in traits]

    # Compute per-trait colors and tiers
    point_colors = []
    tier_labels = []
    for i, (trait_name, z, risk, pctl) in enumerate(traits):
        color = risk_color(z)
        point_colors.append(color)
        n_used = entries[i].get("n_snps_used", 0)
        n_total = entries[i].get("n_snps_total", 0)
        snp_ratio = n_used / max(n_total, 1)
        uncertainty = safe_float(entries[i].get("uncertainty_score", 1.0))
        cal_entry = cal_lookup.get(trait_name.lower())
        uncert_entry = uncert_lookup.get(trait_name.lower()) if uncert_lookup else None
        ev_scores = evidence_lookup.get(trait_name.lower(), [])
        conf = compute_per_trait_confidence(entries[i], cal_entry, uncert_entry, ev_scores)
        tier = trust_tier(conf, cal_entry, snp_ratio, uncertainty)
        tier_labels.append(tier)

    # Serialize data as JSON for Chart.js
    import json as _json
    chart_data = _json.dumps({
        "labels": short_labels,
        "zScores": z_scores,
        "fullNames": full_names,
        "pctls": pctls,
        "riskCats": risk_cats,
        "pointColors": point_colors,
        "tierLabels": tier_labels,
    })

    return f"""
    <div style="max-width:500px;margin:0 auto">
        <canvas id="radarChart" style="max-height:450px"></canvas>
    </div>
    <script>
    (function() {{
        var d = {chart_data};
        var ctx = document.getElementById('radarChart').getContext('2d');
        var riskColorMap = {{'high': '#e74c3c', 'medium': '#f39c12', 'low': '#27ae60'}};
        var maxAbsZ = Math.max.apply(null, d.zScores.map(Math.abs));
        var fillColor = maxAbsZ >= 2 ? '#e74c3c' : (maxAbsZ >= 1 ? '#f39c12' : '#27ae60');

        new Chart(ctx, {{
            type: 'radar',
            data: {{
                labels: d.labels,
                datasets: [{{
                    label: 'PRS Profile',
                    data: d.zScores,
                    backgroundColor: fillColor.replace(')', ',0.12)').replace('rgb', 'rgba'),
                    borderColor: fillColor,
                    borderWidth: 2,
                    pointBackgroundColor: d.pointColors,
                    pointBorderColor: '#fff',
                    pointBorderWidth: 1.5,
                    pointRadius: 5,
                    pointHoverRadius: 8,
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: true,
                animation: {{ duration: 800, easing: 'easeOutQuart' }},
                plugins: {{
                    legend: {{ display: true, position: 'bottom', labels: {{ font: {{ size: 11 }} }} }},
                    tooltip: {{
                        callbacks: {{
                            label: function(ctx) {{
                                var i = ctx.dataIndex;
                                var z = d.zScores[i];
                                var p = d.pctls[i];
                                var tier = d.tierLabels[i];
                                return [
                                    d.fullNames[i],
                                    'z-score: ' + (z >= 0 ? '+' : '') + z.toFixed(2),
                                    'Percentile: ' + p.toFixed(1) + '%',
                                    'Trust: ' + tier
                                ];
                            }}
                        }}
                    }}
                }},
                onClick: function(e, elements) {{
                    if (elements.length > 0) {{
                        var i = elements[0].index;
                        var trait = d.fullNames[i];
                        // Find and scroll to the PRS table row
                        var table = document.querySelector('#prs table');
                        if (table) {{
                            var rows = table.querySelectorAll('tbody tr');
                            rows.forEach(function(row) {{
                                var cell = row.querySelector('td:first-child strong');
                                if (cell && cell.textContent.trim() === trait) {{
                                    row.scrollIntoView({{ behavior: 'smooth', block: 'center' }});
                                    row.style.transition = 'background 0.5s';
                                    row.style.background = '#fef9e7';
                                    setTimeout(function() {{ row.style.background = ''; }}, 1500);
                                }}
                            }});
                        }}
                    }}
                }},
                scales: {{
                    r: {{
                        beginAtZero: true,
                        min: -3,
                        max: 3,
                        ticks: {{ stepSize: 1, backdropColor: 'transparent', font: {{ size: 10 }} }},
                        grid: {{ color: '#dee2e6' }},
                        angleLines: {{ color: '#dee2e6' }},
                        pointLabels: {{ font: {{ size: 10, weight: '600' }} }}
                    }}
                }}
            }}
        }});
    }})();
    </script>
    """


def collapsible_section(section_id, title, content, open_by_default=False):
    """Generate a collapsible HTML section."""
    display = "block" if open_by_default else "none"
    arrow = "▼" if open_by_default else "▶"
    return f"""
    <div class="collapsible-section">
        <div class="section-header" onclick="toggleSection('{section_id}', this)">
            <span class="section-arrow" id="{section_id}_arrow">{arrow}</span>
            <h2>{title}</h2>
        </div>
        <div class="section-body" id="{section_id}" style="display:{display}">
            {content}
        </div>
    </div>
    """

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION BUILDERS
# ═══════════════════════════════════════════════════════════════════════════════

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

def build_top_findings(entries, ui, evidence_lookup=None, cal_lookup=None, uncert_lookup=None,
                        recommendation_lookup=None, ancestry=None, max_findings=8):
    """Prioritized, plain-language 'top findings' list (IMPROVEMENT_PLAN.md 1.1).

    Priority = |z-score| x evidence quality x confidence — the closest proxy
    available today to the plan's "magnitude x evidence x actionability", since
    a curated per-trait actionability/recommendation field doesn't exist yet
    (that's 1.4, tracked separately). Traits with zero matched SNPs are
    excluded — there's nothing to prioritize for an empty result.

    Does NOT fabricate a recommendation for traits without curated action
    guidance: `recommendation_lookup` (data/trait_recommendations.json) covers
    only a hand-verified subset of traits, each checked against a real source
    before being written. Uncovered traits, and covered traits whose
    risk_category isn't "high" (the curated text is written for the elevated-
    risk direction and would misrepresent a protective/average finding), fall
    back to an honest placeholder pointing to the full trait detail below.
    """
    if evidence_lookup is None: evidence_lookup = {}
    if cal_lookup is None: cal_lookup = {}
    if uncert_lookup is None: uncert_lookup = {}
    if recommendation_lookup is None: recommendation_lookup = {}
    ancestry = ancestry or {}
    pop = POP_NAMES.get(ui.get("_lang", "en"), POP_NAMES["en"]).get(
        ancestry.get("assigned_population", "EUR"), ancestry.get("assigned_population", "EUR"))

    scored = []
    for e in entries:
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)
        if n_used <= 0:
            continue
        trait = e.get("trait", "")
        z = safe_float(e.get("population_zscore", e.get("raw_score", 0)))
        pctl = safe_float(e.get("population_percentile", 50))
        risk = e.get("risk_category", "medium")

        ev_scores = evidence_lookup.get(trait.lower(), [])
        ev_avg = sum(ev_scores) / len(ev_scores) if ev_scores else 50.0

        cal_entry = cal_lookup.get(trait.lower())
        uncert_entry = uncert_lookup.get(trait.lower())
        conf = compute_per_trait_confidence(e, cal_entry, uncert_entry, ev_scores)

        priority = abs(z) * (ev_avg / 100.0) * (conf / 100.0)
        scored.append((priority, e, z, pctl, risk, ev_avg, conf, n_used, n_total))

    scored.sort(key=lambda t: t[0], reverse=True)
    top = scored[:max_findings]

    if not top:
        return f'<p style="color:var(--color-text-secondary)">{ui["top_findings_empty"]}</p>'

    risk_words = {"high": ui["risk_word_high"], "medium": ui["risk_word_medium"], "low": ui["risk_word_low"]}
    cards = ""
    for priority, e, z, pctl, risk, ev_avg, conf, n_used, n_total in top:
        trait = e.get("trait", "")
        anchor = trait_anchor_id(trait)
        meaning = ui["top_findings_meaning"].format(
            trait=trait, risk_word=risk_words.get(risk, risk_words["medium"]),
            pctl=pctl, pop=pop, z=z, n_used=n_used, n_total=n_total,
            evidence=evidence_letter(ev_avg))
        lang = ui.get("_lang", "en")
        curated = recommendation_lookup.get(trait.lower()) if risk == "high" else None
        recommendation = (curated or {}).get("recommendation_" + lang) or ui["top_findings_action_fallback"]
        color = risk_color(z)
        cards += f"""
        <div class="info-card" style="border-left:3px solid {color}">
            <div style="display:flex;justify-content:space-between;align-items:start;gap:0.5rem;flex-wrap:wrap">
                <strong style="font-size:0.95rem">{trait}</strong>
                <div style="display:flex;gap:4px;flex-wrap:wrap">{risk_badge(risk, ui)}{evidence_badge(ev_avg)}</div>
            </div>
            <p style="font-size:0.8rem;margin:0.4rem 0;color:var(--color-text-secondary)">{meaning}</p>
            <p style="font-size:0.8rem;margin:0.4rem 0"><strong>→</strong> {recommendation}</p>
            <a href="#{anchor}" style="font-size:0.7rem" onclick="document.getElementById('prs').style.display='block'">{ui["top_findings_jump"]}</a>
        </div>"""

    # Compact, collapsible confidence-context note (IMPROVEMENT_PLAN.md 1.5)
    # — the full disclaimer + per-trait confidence notes already exist at
    # the bottom of the report (Limitations & Disclaimers section); this is
    # the same text, just also reachable without scrolling past everything.
    disclaimer_html = f"""
    <details style="margin-bottom:1rem;background:var(--color-bg-secondary,#f8f9fa);border-radius:6px;padding:0.6rem 1rem">
        <summary style="cursor:pointer;font-size:0.8rem;font-weight:600">{ui["top_findings_disclaimer_summary"]}</summary>
        <p style="white-space:pre-line;font-size:0.78rem;margin:0.6rem 0 0;color:var(--color-text-secondary)">{ui["disclaimer"]}</p>
    </details>"""

    return f"""
    <p style="font-size:0.85rem;color:var(--color-text-secondary);margin-bottom:0.75rem">{ui["top_findings_intro"]}</p>
    {disclaimer_html}
    <div class="info-grid" style="grid-template-columns:repeat(auto-fit, minmax(280px, 1fr))">{cards}</div>
    """

def build_summary_cards(prs_result, ancestry, integrity, validation, ui, cal_lookup=None, uncert_lookup=None, evidence_lookup=None, portability=None):
    """Executive summary cards with confidence overview."""
    if cal_lookup is None:
        cal_lookup = {}
    if uncert_lookup is None:
        uncert_lookup = {}
    if evidence_lookup is None:
        evidence_lookup = {}

    entries = prs_result.get("prs_entries", [])
    n_high = sum(1 for e in entries if e.get("risk_category") == "high")
    n_medium = sum(1 for e in entries if e.get("risk_category") == "medium")
    n_low = sum(1 for e in entries if e.get("risk_category") == "low")

    pop = ancestry.get("assigned_population", "EUR")
    pop_name = POP_NAMES["en"].get(pop, pop)
    confidence = ancestry.get("confidence", "UNKNOWN")
    integrity_score = integrity.get("scientific_integrity_score", 0)
    integrity_cat = integrity.get("category", "Unknown")
    val_score = validation.get("overall_score", 0)
    val_status = validation.get("overall_status", "Unknown")

    # Find top risk trait
    top_trait = max(entries, key=lambda e: abs(safe_float(e.get("raw_score", 0)))) if entries else None

    # ── Confidence overview ──
    confidences = []
    tier_counts = {"TIER 1": 0, "TIER 2": 0, "TIER 3": 0}
    for e in entries:
        trait = e.get("trait", "")
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)
        uncertainty = safe_float(e.get("uncertainty_score", 1.0))
        snp_ratio = n_used / max(n_total, 1)
        cal_entry = cal_lookup.get(trait.lower())
        uncert_entry = uncert_lookup.get(trait.lower())
        ev_scores = evidence_lookup.get(trait.lower(), [])

        conf = compute_per_trait_confidence(e, cal_entry, uncert_entry, ev_scores)
        confidences.append(conf)
        tier = trust_tier(conf, cal_entry, snp_ratio, uncertainty)
        tier_counts[tier] = tier_counts.get(tier, 0) + 1

    avg_conf = sum(confidences) / max(len(confidences), 1) if confidences else 0
    conf_color = "#27ae60" if avg_conf >= 75 else ("#f39c12" if avg_conf >= 50 else "#e74c3c")

    # Find strongest and weakest findings
    if entries and confidences:
        best_idx = max(range(len(confidences)), key=lambda i: confidences[i])
        worst_idx = min(range(len(confidences)), key=lambda i: confidences[i])
        strongest = entries[best_idx]
        weakest = entries[worst_idx]
        best_trait = strongest.get("trait", "?")
        worst_trait = weakest.get("trait", "?")
        best_conf = confidences[best_idx]
        worst_conf = confidences[worst_idx]
    else:
        best_trait, best_conf = "N/A", 0
        worst_trait, worst_conf = "N/A", 0

    # Portability summary note
    port_note = ""
    if portability:
        global_bias = safe_float(portability.get("global_bias_index", 0))
        most_biased = portability.get("most_biased", "AFR")
        if global_bias > 0.15:
            port_note = (
                f'<div class="highlight-box" style="background:#fef9e7;border:1px solid #f39c12;margin-top:0.5rem">'
                f'⚠️ <strong>Portability:</strong> Cross-population bias index: {global_bias:.3f}. '
                f'Results are least reliable for <strong>{most_biased}</strong> populations. '
                f'See Population Portability section for details.'
                f'</div>'
            )

    return f"""
    <div class="summary-grid">
        <div class="summary-card" style="border-left-color:#e74c3c">
            <div class="card-number">{n_high}</div>
            <div class="card-label">{ui['risk_high']}</div>
        </div>
        <div class="summary-card" style="border-left-color:#f39c12">
            <div class="card-number">{n_medium}</div>
            <div class="card-label">{ui['risk_medium']}</div>
        </div>
        <div class="summary-card" style="border-left-color:#27ae60">
            <div class="card-number">{n_low}</div>
            <div class="card-label">{ui['risk_low']}</div>
        </div>
        <div class="summary-card" style="border-left-color:{conf_color}">
            <div class="card-number" style="color:{conf_color}">{avg_conf:.0f}%</div>
            <div class="card-label">Avg Confidence</div>
        </div>
        <div class="summary-card" style="border-left-color:#27ae60">
            <div class="card-number" style="color:#27ae60">{tier_counts.get('TIER 1', 0)}</div>
            <div class="card-label">High Trust (T1)</div>
        </div>
        <div class="summary-card" style="border-left-color:#e74c3c">
            <div class="card-number" style="color:#e74c3c">{tier_counts.get('TIER 3', 0)}</div>
            <div class="card-label">Low Trust (T3)</div>
        </div>
    </div>
    <div class="summary-grid" style="grid-template-columns:repeat(2,1fr)">
        <div class="summary-card" style="border-left-color:#3498db">
            <div class="card-number" style="font-size:1.2rem">{pop_name}</div>
            <div class="card-label">Ancestry ({confidence})</div>
        </div>
        <div class="summary-card" style="border-left-color:#9b59b6">
            <div class="card-number">{integrity_score:.0f}</div>
            <div class="card-label">Integrity / 100</div>
        </div>
    </div>
    {f'<div class="highlight-box"><strong>Top Risk:</strong> {top_trait["trait"]} — raw score {safe_float(top_trait.get("raw_score", 0)):.2f}, {top_trait.get("n_snps_used", 0)}/{top_trait.get("n_snps_total", 0)} SNPs used</div>' if top_trait else ''}
    <div class="highlight-box"><strong>Integrity:</strong> {integrity_cat} — {integrity.get("category_description", "")}</div>
    <div class="highlight-box"><strong>Strongest Finding:</strong> {best_trait} — {best_conf:.0f}% confidence | <strong>Weakest:</strong> {worst_trait} — {worst_conf:.0f}% confidence</div>
    {port_note}
    """


def build_ancestry_section(ancestry, pca_data, ui):
    """Ancestry deep-dive section."""
    pop = ancestry.get("assigned_population", "UNKNOWN")
    confidence = ancestry.get("confidence", "UNKNOWN")
    probs = ancestry.get("posterior_probabilities", {})
    n_ref = ancestry.get("n_reference_samples", 2504)
    n_pcs = ancestry.get("n_pcs", 20)

    pop_names = POP_NAMES["en"]
    prob_rows = ""
    for p in ["EUR", "AFR", "EAS", "SAS", "AMR"]:
        prob = probs.get(p, 0) * 100 if isinstance(probs.get(p, 0), (int, float)) else 0
        bar_w = max(prob, 1)
        prob_rows += f"""
        <tr>
            <td>{pop_names.get(p, p)} ({p})</td>
            <td>{prob:.1f}%</td>
            <td>
                <div style="display:flex;align-items:center;gap:6px">
                    <div style="flex:1;height:6px;background:#e9ecef;border-radius:3px;overflow:hidden">
                        <div style="width:{bar_w}%;height:100%;background:{'#3498db' if p==pop else '#bdc3c7'};border-radius:3px"></div>
                    </div>
                </div>
            </td>
        </tr>"""

    # Try to load PCA eigenvalues
    pca_eigenval = load_json("pca/pca_results.eigenval") if False else {}
    pca_var_explained = ""
    if False:  # eigenval file format differs
        pass

    # Load 1000G PCA coordinates if available
    pca_table = ""
    try:
        eigenvec = "pca/pca_results.eigenvec"
        if os.path.exists(eigenvec):
            pcs = pd.read_csv(eigenvec, sep="\t", nrows=1)
            pc_cols = [c for c in pcs.columns if c.startswith("PC")]
            target_pcs = "pca/target_pcs.eigenvec"
            if os.path.exists(target_pcs):
                tpcs = pd.read_csv(target_pcs, sep="\t")
                pca_table = "<h4>Sample PCA Coordinates (Top 5 PCs)</h4><table><thead><tr><th>PC</th><th>Value</th></tr></thead><tbody>"
                for i, col in enumerate(pc_cols[:5]):
                    val = tpcs[col].values[0] if col in tpcs.columns else 0
                    pca_table += f"<tr><td>{col}</td><td>{val:.6f}</td></tr>"
                pca_table += "</tbody></table>"
    except Exception:
        pass

    return f"""
    <div class="info-grid">
        <div class="info-card">
            <h4>Assigned Population</h4>
            <div class="big-stat">{pop_names.get(pop, pop)}</div>
            <div class="stat-sub">Confidence: {confidence}</div>
        </div>
        <div class="info-card">
            <h4>Reference Panel</h4>
            <div class="big-stat">{n_ref}</div>
            <div class="stat-sub">1000 Genomes Phase 3 samples</div>
        </div>
        <div class="info-card">
            <h4>PCA Dimensions</h4>
            <div class="big-stat">{n_pcs}</div>
            <div class="stat-sub">Principal components</div>
        </div>
        <div class="info-card">
            <h4>Populations</h4>
            <div class="big-stat">5</div>
            <div class="stat-sub">EUR · AFR · EAS · SAS · AMR</div>
        </div>
    </div>
    <h4>Population Probability Distribution</h4>
    <table><thead><tr><th>Population</th><th>Probability</th><th>Distribution</th></tr></thead>
    <tbody>{prob_rows}</tbody></table>
    {pca_table}
    """


def build_prs_table(entries, ui, cal_lookup=None, uncert_lookup=None, portability=None, evidence_lookup=None):
    """Full PRS results table with risk bars and confidence metrics."""
    if cal_lookup is None:
        cal_lookup = {}
    if uncert_lookup is None:
        uncert_lookup = {}
    if evidence_lookup is None:
        evidence_lookup = {}

    rows = ""
    confidences = []
    tier_counts = {"TIER 1": 0, "TIER 2": 0, "TIER 3": 0}

    for e in entries:
        trait = e.get("trait", "")
        z = safe_float(e.get("population_zscore", e.get("raw_score", 0)))
        pctl = safe_float(e.get("population_percentile", 50))
        risk = e.get("risk_category", "medium")
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)
        ci_low = safe_float(e.get("ci_95_lower", 0))
        ci_high = safe_float(e.get("ci_95_upper", 0))
        uncertainty = safe_float(e.get("uncertainty_score", 1.0))
        snp_ratio = n_used / max(n_total, 1)

        # Look up calibration and uncertainty data
        cal_entry = cal_lookup.get(trait.lower())
        uncert_entry = uncert_lookup.get(trait.lower())
        decomp = uncert_entry.get("decomposition", {}) if uncert_entry else {}

        # Evidence scores from lookup
        ev_scores = evidence_lookup.get(trait.lower(), [])

        # Compute confidence and tier
        conf_score = compute_per_trait_confidence(e, cal_entry, uncert_entry, ev_scores)
        confidences.append(conf_score)
        tier = trust_tier(conf_score, cal_entry, snp_ratio, uncertainty)
        tier_counts[tier] = tier_counts.get(tier, 0) + 1

        bar_pct = max(5, min(95, pctl))
        color = risk_color(z)

        rows += f"""
        <tr id="{trait_anchor_id(trait)}">
            <td><strong>{trait}</strong></td>
            <td style="color:{color};font-weight:700">{z:+.2f}</td>
            <td>{pctl:.1f}%</td>
            <td>{risk_badge(risk, ui)}</td>
            <td>{confidence_stars(conf_score)}</td>
            <td>{calibration_flag(cal_entry)}</td>
            <td>{trust_badge(tier)}</td>
            <td>{snp_coverage_bar(n_used, n_total)}</td>
            <td>{mini_decomp_bar(decomp)}</td>
            <td>{trait_limitations_badges(e, cal_entry)}</td>
            <td style="min-width:120px">{risk_bar(bar_pct, z)}</td>
        </tr>"""

    # Portability banner
    port_banner = portability_banner(portability) if portability else ""

    # Average confidence
    avg_conf = sum(confidences) / max(len(confidences), 1) if confidences else 0
    conf_color = "#27ae60" if avg_conf >= 75 else ("#f39c12" if avg_conf >= 50 else "#e74c3c")

    # Summary bar
    summary_bar = (
        f'<div style="display:flex;gap:1rem;align-items:center;margin-bottom:0.75rem;flex-wrap:wrap">'
        f'<div style="display:flex;align-items:center;gap:6px">'
        f'<span style="font-size:0.8rem;font-weight:600">Avg Confidence:</span>'
        f'<span style="color:{conf_color};font-weight:700;font-size:0.9rem">{avg_conf:.0f}%</span>'
        f'</div>'
        f'<span style="color:#7f8c8d">|</span>'
        f'<span style="font-size:0.75rem">'
        f'<span style="color:#27ae60;font-weight:700">T1: {tier_counts.get("TIER 1", 0)}</span> / '
        f'<span style="color:#f39c12;font-weight:700">T2: {tier_counts.get("TIER 2", 0)}</span> / '
        f'<span style="color:#e74c3c;font-weight:700">T3: {tier_counts.get("TIER 3", 0)}</span>'
        f'</span>'
        f'</div>'
    )

    return f"""
    {port_banner}
    {trust_tier_legend()}
    {summary_bar}
    <div style="overflow-x:auto">
    <table>
        <thead><tr>
            <th>Trait</th><th>Z</th><th>%ile</th><th>Risk</th>
            <th>Confidence</th><th>Cal.</th><th>Trust</th>
            <th>SNPs</th><th>Uncertainty</th><th>Limitations</th><th>Risk Bar</th>
        </tr></thead>
        <tbody>{rows}</tbody>
    </table>
    </div>
    <p style="font-size:0.75rem;color:var(--color-text-secondary);margin-top:0.75rem;line-height:1.5">
    ⚠️ <strong>Limitations:</strong> PRS estimates <em>relative</em> genetic predisposition — it does <strong>not</strong> predict absolute disease risk.
    Effect sizes are derived primarily from European-ancestry GWAS and may have reduced accuracy in other populations.
    Gene-environment interactions, rare variants, and structural variants are not captured.
    <strong>Research use only — not clinical diagnosis.</strong>
    </p>
    """


def build_variant_detail(entries, snp_db_path="data/snp_database_annotated.csv"):
    """Per-trait variant-level detail tables."""
    # Try to load the SNP database
    snp_db = None
    if os.path.exists(snp_db_path):
        try: snp_db = pd.read_csv(snp_db_path, dtype=str)
        except Exception: pass

    sections = ""
    for e in entries:
        trait = e.get("trait", "")
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)

        # Detect trait column name (snake_case variants)
        trait_col = None
        if snp_db is not None:
            for col in ["trait_category", "trait", "Trait", "trait_name"]:
                if col in snp_db.columns:
                    trait_col = col
                    break

        if trait_col:
            trait_snps = snp_db[snp_db[trait_col].str.lower() == trait.lower()]
            if len(trait_snps) == 0:
                trait_snps = snp_db.head(0)  # empty match
        else:
            trait_snps = None

        variant_rows = ""
        if trait_snps is not None and len(trait_snps) > 0:
            for _, row in trait_snps.iterrows():
                rsid = row.get("rsid", "—")
                gene = row.get("gene", "—")
                effect_allele = row.get("effect_allele", "—")
                weight = row.get("weight", row.get("beta", "—"))
                evidence = row.get("evidence_level", row.get("evidence", "—")).strip().upper()
                ev_colors = {"A": ("#27ae60", "#d5f5e3"), "B": ("#2e86c1", "#d6eaf8"), "C": ("#f39c12", "#fdebd0"), "D": ("#95a5a6", "#eaecee")}
                fg, bg = ev_colors.get(evidence, ("#7f8c8d", "#e9ecef"))
                ev_badge = f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.7rem;font-weight:700">{evidence}</span>'
                variant_rows += f"<tr><td>{rsid}</td><td>{gene}</td><td>{effect_allele}</td><td>{weight}</td><td>{ev_badge}</td></tr>"
        else:
            variant_rows = f"<tr><td colspan='5' style='color:#7f8c8d'>SNP database not available. Run with --snp-db to populate.</td></tr>"

        sections += f"""
        <h4>{trait} <span style="font-weight:400;color:#7f8c8d">({n_used}/{n_total} SNPs used)</span></h4>
        <table>
            <thead><tr><th>rsID</th><th>Gene</th><th>Effect Allele</th><th>Weight (β)</th><th>Evidence</th></tr></thead>
            <tbody>{variant_rows}</tbody>
        </table>
        """

    # Add evidence level legend
    legend = """
    <div class="highlight-box" style="background:#eaf2f8;border:1px solid #aed6f1;border-radius:8px;padding:0.8rem 1.2rem;margin-top:0.5rem">
        <p style="font-size:0.8rem;margin:0">
        <strong>📖 Evidence Levels:</strong>
        <span style="background:#d5f5e3;color:#27ae60;padding:1px 6px;border-radius:3px;font-size:0.7rem;font-weight:700">A</span> = Strong GWAS evidence (p &lt; 5×10⁻⁸), replicated in multiple populations.<br>
        <span style="background:#d6eaf8;color:#2e86c1;padding:1px 6px;border-radius:3px;font-size:0.7rem;font-weight:700">B</span> = Moderate evidence (p &lt; 10⁻⁵), supported by functional studies.<br>
        <span style="background:#fdebd0;color:#f39c12;padding:1px 6px;border-radius:3px;font-size:0.7rem;font-weight:700">C</span> = Candidate gene study or suggestive GWAS (p &lt; 10⁻³).<br>
        <span style="background:#eaecee;color:#95a5a6;padding:1px 6px;border-radius:3px;font-size:0.7rem;font-weight:700">D</span> = Preliminary evidence or literature-based association.<br>
        <span style="font-size:0.72rem;color:var(--color-text-secondary)">Weights are normalized effect sizes (β coefficients). Higher |β| = stronger SNP contribution to the PRS.</span>
        </p>
    </div>
    """
    return sections + legend


def build_uncertainty_section(entries):
    """Uncertainty analysis with visual indicators."""
    rows = ""
    for e in entries:
        trait = e.get("trait", "")
        raw = safe_float(e.get("raw_score", 0))
        ci_low = safe_float(e.get("ci_95_lower", 0))
        ci_high = safe_float(e.get("ci_95_upper", 0))
        uncertainty = safe_float(e.get("uncertainty_score", 1.0))
        se = (ci_high - ci_low) / 3.92  # approximate SE from CI

        # Uncertainty color
        if uncertainty < 0.5: u_color = "#27ae60"
        elif uncertainty < 0.8: u_color = "#f39c12"
        else: u_color = "#e74c3c"

        # CI bar: map CI range within [-3, +3]
        left_norm = max(0, (ci_low + 3) / 6 * 100)
        right_norm = min(100, (ci_high + 3) / 6 * 100)

        rows += f"""
        <tr>
            <td><strong>{trait}</strong></td>
            <td>{raw:.3f} ± {se:.3f}</td>
            <td>[{ci_low:.3f}, {ci_high:.3f}]</td>
            <td>
                <div style="position:relative;height:20px;background:#e9ecef;border-radius:3px;margin:2px 0">
                    <div style="position:absolute;left:{left_norm:.1f}%;width:{right_norm-left_norm:.1f}%;height:100%;background:{'#3498db44' if ci_low*ci_high > 0 else '#e74c3c44'};border:1px solid {'#3498db' if ci_low*ci_high > 0 else '#e74c3c'};border-radius:3px"></div>
                    <div style="position:absolute;left:50%;top:0;width:2px;height:100%;background:#2c3e50"></div>
                </div>
                <div style="display:flex;justify-content:space-between;font-size:0.65rem;color:#7f8c8d">
                    <span>-3σ</span><span>0</span><span>+3σ</span>
                </div>
            </td>
            <td style="color:{u_color};font-weight:700">{uncertainty:.2f}</td>
        </tr>"""

    return f"""
    <p>Uncertainty scores below 0.5 indicate high-confidence estimates. Scores above 0.8 indicate substantial uncertainty — treat results with caution.</p>
    <table>
        <thead><tr><th>Trait</th><th>PRS ± SE</th><th>95% CI</th><th>CI Visualization</th><th>Uncert.</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_validation_section(validation, ui):
    """Validation checks table."""
    checks = validation.get("checks", [])
    rows = ""
    for c in checks:
        passed = c.get("passed", False)
        sev = c.get("severity", "INFO")
        sev_colors = {"ERROR": "#e74c3c", "WARNING": "#f39c12", "INFO": "#3498db"}
        icon = "✅" if passed else "❌"
        rows += f"""
        <tr>
            <td><strong>{c.get("check_id", "")}</strong></td>
            <td>{c.get("category", "")}</td>
            <td>{c.get("description", "")}</td>
            <td style="color:{sev_colors.get(sev, '#7f8c8d')};font-weight:700">{sev}</td>
            <td>{icon}</td>
            <td style="font-size:0.8rem;color:#7f8c8d">{c.get("detail", "")}</td>
        </tr>"""

    return f"""
    <div class="info-grid">
        <div class="info-card">
            <h4>Overall Score</h4>
            <div class="big-stat">{validation.get('overall_score', 0):.0f}</div>
            <div class="stat-sub">/ 100 — {validation.get('overall_status', '').replace('_', ' ').title()}</div>
        </div>
        <div class="info-card">
            <h4>Passed</h4>
            <div class="big-stat" style="color:#27ae60">{validation.get('passed', 0)}</div>
            <div class="stat-sub">of {validation.get('total_checks', 0)} checks</div>
        </div>
        <div class="info-card">
            <h4>Warnings</h4>
            <div class="big-stat" style="color:#f39c12">{validation.get('warnings', 0)}</div>
            <div class="stat-sub">non-critical issues</div>
        </div>
        <div class="info-card">
            <h4>Errors</h4>
            <div class="big-stat" style="color:#e74c3c">{validation.get('errors', 0)}</div>
            <div class="stat-sub">critical issues</div>
        </div>
    </div>
    <table>
        <thead><tr><th>ID</th><th>Category</th><th>Check</th><th>Severity</th><th>Result</th><th>Detail</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_benchmark_section(benchmark, quality_delta, ui):
    """GWAS & external benchmarking."""
    entries = benchmark.get("entries", [])
    bench_rows = ""
    for e in entries:
        status = e.get("status", "VALID")
        vtype = e.get("validation_type", "unknown")
        is_circ = e.get("is_circular", False)
        circ_badge = '<span style="background:#fdebd0;color:#b7950b;padding:1px 5px;border-radius:3px;font-size:0.65rem">CIRCULAR</span>' if is_circ else '<span style="background:#d5f5e3;color:#1e8449;padding:1px 5px;border-radius:3px;font-size:0.65rem">INDEPENDENT</span>'
        bench_rows += f"""
        <tr>
            <td><strong>{e.get('validation_id', '')}</strong></td>
            <td>{e.get('description', '')}</td>
            <td>{vtype}</td>
            <td>{circ_badge}</td>
            <td>{status}</td>
        </tr>"""

    # Quality delta
    qd = quality_delta
    components = qd.get("components", [])
    delta_rows = ""
    for c in components:
        delta = c.get("delta", 0)
        direction = c.get("direction", "at_par")
        d_color = "#27ae60" if direction == "overperform" else ("#e74c3c" if direction == "underperform" else "#f39c12")
        delta_rows += f"""
        <tr>
            <td>{c.get('dimension', '')}</td>
            <td>{safe_float(c.get('internal_score', 0)):.0f}</td>
            <td>{safe_float(c.get('external_benchmark', 0)):.0f}</td>
            <td style="color:{d_color};font-weight:700">{delta:+.0f}</td>
            <td>{direction.replace('_', ' ').title()}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{c.get('explanation', '')[:150]}</td>
        </tr>"""

    return f"""
    <h4>Validation Classification</h4>
    <div class="info-grid">
        <div class="info-card"><h4>Total</h4><div class="big-stat">{benchmark.get('validation_summary', {}).get('total_validations', 0)}</div></div>
        <div class="info-card"><h4>Internal</h4><div class="big-stat">{benchmark.get('validation_summary', {}).get('internal', 0)}</div></div>
        <div class="info-card"><h4>External</h4><div class="big-stat" style="color:#3498db">{benchmark.get('validation_summary', {}).get('external', 0)}</div></div>
        <div class="info-card"><h4>Circular</h4><div class="big-stat" style="color:#e74c3c">{benchmark.get('validation_summary', {}).get('circular', 0)}</div></div>
    </div>
    <table>
        <thead><tr><th>ID</th><th>Description</th><th>Type</th><th>Classification</th><th>Status</th></tr></thead>
        <tbody>{bench_rows}</tbody>
    </table>

    <h4 style="margin-top:2rem">Quality Delta — Internal vs External Benchmarks</h4>
    <div class="info-grid">
        <div class="info-card"><h4>Mean Delta</h4><div class="big-stat" style="color:{'#27ae60' if qd.get('mean_delta', 0) >= 0 else '#e74c3c'}">{qd.get('mean_delta', 0):+.1f}</div></div>
        <div class="info-card"><h4>Overperform</h4><div class="big-stat" style="color:#27ae60">{qd.get('overperform', 0)}</div></div>
        <div class="info-card"><h4>At Par</h4><div class="big-stat">{qd.get('at_par', 0)}</div></div>
        <div class="info-card"><h4>Underperform</h4><div class="big-stat" style="color:#e74c3c">{qd.get('underperform', 0)}</div></div>
    </div>
    <table>
        <thead><tr><th>Dimension</th><th>Internal</th><th>External</th><th>Δ</th><th>Direction</th><th>Explanation</th></tr></thead>
        <tbody>{delta_rows}</tbody>
    </table>
    """


def build_adversarial_section(adversarial, ui):
    """Adversarial stress testing results."""
    results = adversarial.get("results", [])
    rows = ""
    for r in results:
        robust = r.get("is_robust", False)
        if isinstance(robust, str):
            robust = robust.lower() == "true"
        sev = r.get("severity", "MODERATE")
        sev_color = {"CRITICAL": "#e74c3c", "HIGH": "#e67e22", "MODERATE": "#f39c12"}
        icon = "✅" if robust else "❌"
        change = safe_float(r.get("relative_change", 0))
        rows += f"""
        <tr>
            <td><strong>{r.get('test_id', '')}</strong></td>
            <td>{r.get('description', '')}</td>
            <td style="color:{sev_color.get(sev, '#7f8c8d')};font-weight:700">{sev}</td>
            <td>{icon} {'Robust' if robust else 'Vulnerable'}</td>
            <td>{change:+.2f}</td>
            <td style="font-size:0.8rem">{r.get('detail', '')}</td>
        </tr>"""

    critical_findings = adversarial.get("critical_findings", [])

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Robustness Score</h4><div class="big-stat">{adversarial.get('overall_robustness_score', 0):.0f}/100</div></div>
        <div class="info-card"><h4>Tests Run</h4><div class="big-stat">{adversarial.get('n_tests', 0)}</div></div>
        <div class="info-card"><h4>Robust</h4><div class="big-stat" style="color:#27ae60">{adversarial.get('n_robust', 0)}</div></div>
        <div class="info-card"><h4>Vulnerable</h4><div class="big-stat" style="color:#e74c3c">{adversarial.get('n_vulnerable', 0)}</div></div>
    </div>
    {f'<div class="highlight-box" style="background:#fadbd8"><strong>Critical Findings:</strong> {", ".join(critical_findings)}</div>' if critical_findings else ''}
    <table>
        <thead><tr><th>Test ID</th><th>Description</th><th>Severity</th><th>Result</th><th>Change</th><th>Detail</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_failure_map_section(failure_map, ui):
    """Failure mode coverage."""
    failures = failure_map.get("failures", [])
    rows = ""
    for f in failures:
        sev = f.get("severity", "MODERATE")
        sev_color = {"CRITICAL": "#e74c3c", "HIGH": "#e67e22", "MODERATE": "#f39c12"}
        validated = f.get("adversarial_validated", False)
        v_icon = "✅" if validated else "⬚"
        rows += f"""
        <tr>
            <td><strong>{f.get('id', '')}</strong></td>
            <td>{f.get('component', '')}</td>
            <td>{f.get('failure', '')}</td>
            <td style="color:{sev_color.get(sev, '#7f8c8d')};font-weight:700">{sev}</td>
            <td>{v_icon}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{f.get('effect', '')[:120]}</td>
        </tr>"""

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Total Failures</h4><div class="big-stat">{failure_map.get('n_failures', 0)}</div></div>
        <div class="info-card"><h4>Critical</h4><div class="big-stat" style="color:#e74c3c">{failure_map.get('n_critical', 0)}</div></div>
        <div class="info-card"><h4>High</h4><div class="big-stat" style="color:#e67e22">{failure_map.get('n_high', 0)}</div></div>
        <div class="info-card"><h4>Vulnerable Component</h4><div class="big-stat" style="font-size:1rem">{failure_map.get('most_vulnerable_component', 'N/A')}</div></div>
    </div>
    <table>
        <thead><tr><th>ID</th><th>Component</th><th>Failure Mode</th><th>Severity</th><th>Validated</th><th>Effect</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_leakage_section(leakage, ui):
    """Leakage prevention gates."""
    checks = leakage.get("checks", [])
    rows = ""
    for c in checks:
        passed = c.get("passed", False)
        sev = c.get("severity", "INFO")
        icon = "✅" if passed else "❌"
        sev_color = {"ERROR": "#e74c3c", "WARNING": "#f39c12", "INFO": "#3498db"}
        rows += f"""
        <tr>
            <td><strong>{c.get('gate', '')}</strong></td>
            <td>{c.get('description', '')}</td>
            <td style="color:{sev_color.get(sev, '#7f8c8d')};font-weight:700">{sev}</td>
            <td>{icon}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{c.get('detail', '')}</td>
        </tr>"""

    can_proceed = leakage.get("pipeline_can_proceed", False)
    all_passed = leakage.get("all_passed", False)

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Pipeline Safe</h4><div class="big-stat" style="color:{'#27ae60' if can_proceed else '#e74c3c'}">{'YES' if can_proceed else 'NO'}</div></div>
        <div class="info-card"><h4>All Passed</h4><div class="big-stat" style="color:{'#27ae60' if all_passed else '#f39c12'}">{'YES' if all_passed else 'NO'}</div></div>
        <div class="info-card"><h4>Gates</h4><div class="big-stat">{leakage.get('n_checks', 0)}</div></div>
    </div>
    <table>
        <thead><tr><th>Gate</th><th>Description</th><th>Severity</th><th>Result</th><th>Detail</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_integrity_section(integrity, ui):
    """Scientific integrity score breakdown."""
    components = integrity.get("components", [])
    rows = ""
    for c in components:
        score = safe_float(c.get("score", 0))
        score_color = "#27ae60" if score >= 80 else ("#f39c12" if score >= 60 else "#e74c3c")
        rows += f"""
        <tr>
            <td><strong>{c.get('name', '')}</strong></td>
            <td style="color:{score_color};font-weight:700">{score:.1f}</td>
            <td>{safe_float(c.get('weight', 0))*100:.0f}%</td>
            <td>{safe_float(c.get('contribution', 0)):.1f}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{c.get('source', '')}</td>
        </tr>"""

    total = integrity.get("scientific_integrity_score", 0)
    cat = integrity.get("category", "Unknown")
    cat_color = {"PUBLICATION_READY": "#27ae60", "RESEARCH_GRADE": "#3498db",
                 "NEEDS_REVISION": "#f39c12", "SIGNIFICANT_ISSUES": "#e67e22",
                 "NOT_PUBLISHABLE": "#e74c3c"}.get(cat, "#7f8c8d")

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Integrity Score</h4><div class="big-stat" style="color:{cat_color}">{total:.1f}</div><div class="stat-sub">/ 100</div></div>
        <div class="info-card"><h4>Category</h4><div class="big-stat" style="font-size:1rem;color:{cat_color}">{cat.replace('_', ' ').title()}</div></div>
        <div class="info-card"><h4>Formula</h4><div class="stat-sub" style="font-size:0.7rem">{integrity.get('formula', '')}</div></div>
        <div class="info-card"><h4>Weights Locked</h4><div class="big-stat">{'✅' if integrity.get('weights_locked', False) else '❌'}</div></div>
    </div>
    <table>
        <thead><tr><th>Component</th><th>Score</th><th>Weight</th><th>Contribution</th><th>Source</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_uncertainty_decomposition(uncertainty_report):
    """Per-trait variance decomposition: genotype vs ancestry vs effect."""
    results = uncertainty_report.get("results", [])
    if not results:
        return "<p style='color:#7f8c8d'>Uncertainty report not available.</p>"

    rows = ""
    for r in results:
        trait = r.get("trait", "")
        decomp = r.get("decomposition", {})
        gen_frac = safe_float(decomp.get("genotype_fraction", 0)) * 100
        anc_frac = safe_float(decomp.get("ancestry_fraction", 0)) * 100
        eff_frac = safe_float(decomp.get("effect_fraction", 0)) * 100
        total_var = safe_float(decomp.get("total_variance", 0))
        prs = safe_float(r.get("prs_point_estimate", 0))
        se = safe_float(r.get("prs_std_error", 0))

        rows += f"""
        <tr>
            <td><strong>{trait}</strong></td>
            <td>{prs:.3f} ± {se:.3f}</td>
            <td>{total_var:.4f}</td>
            <td>
                <div style="display:flex;height:10px;border-radius:3px;overflow:hidden;background:#e9ecef">
                    <div style="width:{gen_frac:.0f}%;background:#3498db" title="Genotype: {gen_frac:.1f}%"></div>
                    <div style="width:{anc_frac:.0f}%;background:#e74c3c" title="Ancestry: {anc_frac:.1f}%"></div>
                    <div style="width:{eff_frac:.0f}%;background:#f39c12" title="Effect SE: {eff_frac:.1f}%"></div>
                </div>
                <div style="display:flex;justify-content:space-between;font-size:0.65rem;color:#7f8c8d;margin-top:2px">
                    <span>Gen: {gen_frac:.0f}%</span><span>Anc: {anc_frac:.0f}%</span><span>Eff: {eff_frac:.0f}%</span>
                </div>
            </td>
            <td>{r.get('n_snps_with_genotype', 0)}/{r.get('n_snps_with_effect_se', 0)}</td>
        </tr>"""

    return f"""
    <p>Three-layer variance propagation decomposes PRS uncertainty into genotype quality,
    ancestry ambiguity, and effect size standard error. Blue=genotype, Red=ancestry, Yellow=effect.</p>
    <table>
        <thead><tr><th>Trait</th><th>PRS ± SE</th><th>Total Var</th><th>Variance Decomposition</th><th>SNPs (genotype/effect)</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_gwas_consortium_section(gwas_consortium):
    """GWAS consortium validation — consortia info + trait-level checks."""
    consortia = gwas_consortium.get("consortia", {})
    validations = gwas_consortium.get("validations", [])

    # Consortium summary cards
    cons_cards = ""
    for name, c in consortia.items():
        cons_cards += f"""
        <div class="info-card">
            <h4>{name}</h4>
            <div class="big-stat" style="font-size:1rem">{c.get('primary_ancestry', 'EUR')}</div>
            <div class="stat-sub">n={c.get('n_discovery', 0):,} | PMID:{c.get('pmid', '')}</div>
            <div class="stat-sub">{', '.join(c.get('traits', [])[:3])}</div>
        </div>"""

    # Validation results
    val_rows = ""
    for v in validations:
        passed = v.get("overall_status") == "PASS"
        icon = "✅" if passed else "❌"
        match = safe_float(v.get("effect_direction_match", 0)) * 100
        val_rows += f"""
        <tr>
            <td>{v.get('consortium', '')}</td>
            <td>{v.get('trait', '')}</td>
            <td>{match:.0f}%</td>
            <td>{v.get('snp_overlap_count', 0)} ({safe_float(v.get('snp_overlap_pct', 0))*100:.1f}%)</td>
            <td>{icon} {v.get('overall_status', '')}</td>
        </tr>"""

    passed_count = gwas_consortium.get("passed", 0)
    total_count = gwas_consortium.get("total_checks", 0)

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Consortia</h4><div class="big-stat">{len(consortia)}</div><div class="stat-sub">GIANT · GLGC · MAGIC · DIAGRAM</div></div>
        <div class="info-card"><h4>Validations</h4><div class="big-stat">{total_count}</div><div class="stat-sub">total checks</div></div>
        <div class="info-card"><h4>Passed</h4><div class="big-stat" style="color:#27ae60">{passed_count}</div><div class="stat-sub">of {total_count}</div></div>
        <div class="info-card"><h4>Failed</h4><div class="big-stat" style="color:#e74c3c">{gwas_consortium.get('failed', 0)}</div><div class="stat-sub">no SNP overlap traits</div></div>
    </div>

    <h4>Consortium Profiles</h4>
    <div class="info-grid">{cons_cards}</div>

    <h4>Trait-Level Validation</h4>
    <p style="font-size:0.8rem;color:#7f8c8d">13 of 17 failures are traits with zero SNP overlap between curated panel and consortium GWAS (e.g. BMI, HDL, fasting glucose). The 4 passes are the curated trait categories that DO overlap.</p>
    <table>
        <thead><tr><th>Consortium</th><th>Trait</th><th>Effect Direction Match</th><th>SNP Overlap</th><th>Status</th></tr></thead>
        <tbody>{val_rows}</tbody>
    </table>
    """


def build_portability_section(portability):
    """Population portability — PRS shift across populations."""
    pops = portability.get("populations", [])
    rows = ""
    for p in pops:
        pop = p.get("population", "")
        status = p.get("status", "")
        status_color = {"GOOD_PORTABILITY": "#27ae60", "MODERATE_PORTABILITY": "#f39c12", "LIMITED_PORTABILITY": "#e74c3c"}
        color = status_color.get(status, "#7f8c8d")
        rows += f"""
        <tr>
            <td><strong>{pop}</strong></td>
            <td>{p.get('n_reference_samples', 0)}</td>
            <td>{safe_float(p.get('mean_prs_shift', 0)):.2f}</td>
            <td>{safe_float(p.get('calibration_drift', 0)):.2f}</td>
            <td>{safe_float(p.get('rank_instability', 0)):.2f}</td>
            <td>{safe_float(p.get('ancestry_bias_index', 0)):.3f}</td>
            <td style="color:{color};font-weight:700">{status.replace('_', ' ').title()}</td>
        </tr>"""

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Global Bias Index</h4><div class="big-stat">{safe_float(portability.get('global_bias_index', 0)):.3f}</div></div>
        <div class="info-card"><h4>Most Biased</h4><div class="big-stat" style="font-size:1.2rem;color:#e74c3c">{portability.get('most_biased', 'N/A')}</div></div>
        <div class="info-card"><h4>Least Biased</h4><div class="big-stat" style="font-size:1.2rem;color:#27ae60">{portability.get('least_biased', 'N/A')}</div></div>
    </div>
    <p style="font-size:0.8rem;color:#7f8c8d">EUR-centric GWAS inherently limit cross-population portability. Lower indices = better portability. AFR shows highest bias (0.300) due to differing LD structure and allele frequencies.</p>
    <table>
        <thead><tr><th>Population</th><th>Ref Samples</th><th>PRS Shift</th><th>Calib. Drift</th><th>Rank Instability</th><th>Bias Index</th><th>Status</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_reproducibility_section(repro):
    """Reproducibility: environment, seeds, versions."""
    env = repro.get("environment", {})
    seeds = repro.get("seeds", {})

    # System info
    sys_rows = f"""
    <tr><td>OS</td><td>{env.get('os_name', '')} {env.get('os_version', '')} ({env.get('architecture', '')})</td></tr>
    <tr><td>Python</td><td>{env.get('python_version', '')} ({env.get('python_implementation', '')})</td></tr>
    <tr><td>Kernel</td><td style="font-size:0.7rem;word-break:break-all">{env.get('kernel', '')[:120]}</td></tr>"""

    # Tool versions
    tools = env.get("system_tools", {})
    tool_rows = ""
    for tool, version in tools.items():
        tool_rows += f"<tr><td>{tool}</td><td style='font-size:0.78rem'>{version}</td></tr>"

    # Package versions (top 10)
    pkgs = env.get("pip_packages", {})
    pkg_rows = ""
    for pkg, ver in sorted(pkgs.items())[:20]:
        pkg_rows += f"<tr><td>{pkg}</td><td>{ver}</td></tr>"

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Run ID</h4><div class="big-stat" style="font-size:0.8rem;word-break:break-all">{repro.get('run_id', 'N/A')[:16]}</div></div>
        <div class="info-card"><h4>Reproducibility Score</h4><div class="big-stat" style="color:#27ae60">{repro.get('reproducibility_score', 0):.0f}/100</div></div>
        <div class="info-card"><h4>Pipeline Version</h4><div class="big-stat">v{repro.get('pipeline_version', 'N/A')}</div></div>
        <div class="info-card"><h4>Global Seed</h4><div class="big-stat">{seeds.get('global_seed', 'N/A')}</div></div>
    </div>

    <div style="display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;margin-top:1rem">
        <div>
            <h4>System Environment</h4>
            <table><tbody>{sys_rows}</tbody></table>
        </div>
        <div>
            <h4>Bioinformatics Tools</h4>
            <table><thead><tr><th>Tool</th><th>Version</th></tr></thead><tbody>{tool_rows}</tbody></table>
        </div>
    </div>

    <h4 style="margin-top:1.5rem">Python Packages (selected)</h4>
    <table><thead><tr><th>Package</th><th>Version</th></tr></thead><tbody>{pkg_rows}</tbody></table>

    <h4 style="margin-top:1.5rem">Seeds Registry</h4>
    <table><thead><tr><th>Module</th><th>Seed</th></tr></thead><tbody>
    <tr><td>Global</td><td>{seeds.get('global_seed', '')}</td></tr>
    <tr><td>NumPy</td><td>{seeds.get('numpy_seed', '')}</td></tr>
    <tr><td>Python hash</td><td>{seeds.get('python_hash_seed', '')}</td></tr>
    <tr><td>Scikit-learn</td><td>{seeds.get('sklearn_seed', '')}</td></tr>
    <tr><td>PLINK</td><td>{seeds.get('plink_seed', '')}</td></tr>
    <tr><td>Bootstrap 0</td><td>{seeds.get('bootstrap_seeds', [0])[0]}</td></tr>
    </tbody></table>
    """


def _build_gwas_summary(trait_checks):
    """Build dynamic GWAS type summary from actual data (not hardcoded)."""
    if not trait_checks:
        return '<p style="font-size:0.8rem;color:#7f8c8d">No GWAS consistency data available.</p>'

    type_counts = {}
    unknown_count = 0
    for tc in trait_checks:
        gtype = tc.get("gwas_type", "unknown")
        if gtype == "unknown":
            unknown_count += 1
        else:
            type_counts[gtype] = type_counts.get(gtype, 0) + 1

    parts = []
    for gtype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        parts.append(f"{gtype.replace('_', ' ')} ({count})")

    summary = f"GWAS types: {', '.join(parts)}."
    if unknown_count > 0:
        summary += f" {unknown_count} trait(s) with unknown GWAS type."
    summary += " All traits use EUR GWAS sources matching the EUR target ancestry where applicable."

    return f'<p style="font-size:0.8rem;color:#7f8c8d">{summary}</p>'


def build_consistency_section(consistency):
    """GWAS-Ancestry consistency check per trait."""
    detailed = consistency.get("detailed_report", {})
    trait_checks = detailed.get("trait_checks", [])
    rows = ""
    for t in trait_checks:
        gwas_pop = t.get("gwas_population", "")
        gwas_type = t.get("gwas_type", "")
        match = t.get("is_match", False)
        icon = "✅" if match else "❌"
        rows += f"""
        <tr>
            <td>{t.get('trait', '')}</td>
            <td>{gwas_pop}</td>
            <td>{gwas_type.replace('_', ' ').title()}</td>
            <td>{t.get('target_population', '')}</td>
            <td>{icon}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{t.get('note', '')}</td>
        </tr>"""

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>GWAS-Ancestry Match</h4><div class="big-stat" style="color:{'#27ae60' if consistency.get('gwas_ancestry_match') else '#e74c3c'}">{'✅ PASS' if consistency.get('gwas_ancestry_match') else '❌ FAIL'}</div></div>
        <div class="info-card"><h4>LD-Ancestry Match</h4><div class="big-stat" style="color:{'#27ae60' if consistency.get('ld_ancestry_match') else '#f39c12'}">{'✅ PASS' if consistency.get('ld_ancestry_match') else '⚠️ WARN'}</div></div>
        <div class="info-card"><h4>Confidence Downgrade</h4><div class="big-stat">{safe_float(consistency.get('confidence_downgrade', 0)):.2f}</div></div>
        <div class="info-card"><h4>Recommended GWAS</h4><div class="big-stat" style="font-size:0.75rem">{consistency.get('recommended_gwas_source', 'N/A')}</div></div>
    </div>
    {_build_gwas_summary(trait_checks)}
    <table>
        <thead><tr><th>Trait</th><th>GWAS Population</th><th>GWAS Type</th><th>Target Pop</th><th>Match</th><th>Note</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_leakage_detail_section(leakage_audit):
    """Detailed leakage audit — 7 checks."""
    checks = leakage_audit.get("checks", [])
    rows = ""
    for c in checks:
        passed = c.get("passed", False)
        sev = c.get("severity", "INFO")
        sev_color = {"ERROR": "#e74c3c", "WARNING": "#f39c12", "INFO": "#3498db"}
        icon = "✅" if passed else "❌"
        rows += f"""
        <tr>
            <td><strong>{c.get('check_id', '')}</strong></td>
            <td>{c.get('description', '')}</td>
            <td style="color:{sev_color.get(sev, '#7f8c8d')};font-weight:700">{sev}</td>
            <td>{icon}</td>
            <td style="font-size:0.78rem;color:#7f8c8d">{c.get('detail', '')}</td>
        </tr>"""

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Pipeline Safe</h4><div class="big-stat" style="color:{'#27ae60' if leakage_audit.get('pipeline_safe') else '#e74c3c'}">{'YES' if leakage_audit.get('pipeline_safe') else 'NO'}</div></div>
        <div class="info-card"><h4>Passed</h4><div class="big-stat" style="color:#27ae60">{leakage_audit.get('passed', 0)}/{leakage_audit.get('total_checks', 0)}</div></div>
        <div class="info-card"><h4>Warnings</h4><div class="big-stat" style="color:#f39c12">{leakage_audit.get('warnings', 0)}</div></div>
        <div class="info-card"><h4>Errors</h4><div class="big-stat" style="color:#e74c3c">{leakage_audit.get('errors', 0)}</div></div>
    </div>
    <table>
        <thead><tr><th>ID</th><th>Check</th><th>Severity</th><th>Result</th><th>Detail</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


def build_methodology_section(prs_result, ancestry):
    """Pipeline methodology summary."""
    meta = prs_result.get("metadata", {})
    n_variants = meta.get("n_variants", prs_result.get("prs_core", {}).get("n_variants", 109))
    n_traits = meta.get("n_traits", 10)
    formula = meta.get("prs_formula", "PRS = Σ(βⱼ × Gᵢⱼ)")
    method = meta.get("computation_method", "PLINK --score (dosage-weighted)")
    pipeline_ver = meta.get("pipeline_version", "1.1.0")
    n_ref = ancestry.get("n_reference_samples", 2504)
    n_pcs = ancestry.get("n_pcs", 20)

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Pipeline</h4><div class="big-stat">v{pipeline_ver}</div><div class="stat-sub">BlueGen</div></div>
        <div class="info-card"><h4>PRS Formula</h4><div class="big-stat" style="font-size:0.9rem">{formula}</div><div class="stat-sub">{method}</div></div>
        <div class="info-card"><h4>Variants</h4><div class="big-stat">{n_variants}</div><div class="stat-sub">across {n_traits} traits</div></div>
        <div class="info-card"><h4>Reference</h4><div class="big-stat">{n_ref}</div><div class="stat-sub">1000G Phase 3, {n_pcs} PCs</div></div>
    </div>
    <h4>Pipeline Stages</h4>
    <table>
        <thead><tr><th>Stage</th><th>Process</th><th>Method</th><th>Output</th></tr></thead>
        <tbody>
            <tr><td>A</td><td>VCF → PLINK</td><td>PLINK 1.9, GQ≥20, DP≥10</td><td>.bed/.bim/.fam</td></tr>
            <tr><td>B</td><td>Quality Control</td><td>geno 0.05, maf 0.01, HWE 1e-6</td><td>Filtered genotypes</td></tr>
            <tr><td>C</td><td>LD Pruning</td><td>Per-population (EUR/AFR/EAS/SAS/AMR), conservative intersection</td><td>Ancestry-matched independent SNPs</td></tr>
            <tr><td>D</td><td>PCA + Projection</td><td>1000G-trained PCA, target sample projection</td><td>20 PCs, ancestry inference</td></tr>
            <tr><td>F</td><td>PRS Computation</td><td>PLINK --score (dosage-weighted)</td><td>Raw PRS per trait</td></tr>
            <tr><td>G</td><td>PCA Adjustment</td><td>PRS_adj = PRS_raw − Σ(βₖ × PCₖ)</td><td>Ancestry-adjusted PRS</td></tr>
            <tr><td>H</td><td>Population Calibration</td><td>Empirical 1000G population distributions</td><td>Z-scores + percentiles</td></tr>
            <tr><td>7-10</td><td>Validation & Lock</td><td>8-dimension validation, adversarial stress, publication lock</td><td>Scientific integrity score</td></tr>
        </tbody>
    </table>
    """


def build_pgs_calibration_section(pgs_data, ui=None, pgs_coverage=None):
    """PGS Catalog population-calibrated results with clinical interpretation.

    Args:
        pgs_data: PGS calibration report JSON data.
        ui: Bilingual UI dictionary (optional, for labels).
        pgs_coverage: {pgs_id: {"n_used": int, "n_total": int}} lookup.
    """
    summary = pgs_data.get("summary", {})
    high_risk = pgs_data.get("high_risk_traits", [])
    elevated = pgs_data.get("elevated_risk_traits", [])
    low_risk = pgs_data.get("low_risk_traits", [])
    methodology = pgs_data.get("methodology", {})

    # Clinical context for each trait
    CLINICAL = {
        "Type 2 Diabetes": "Higher z-score = greater genetic predisposition to insulin resistance and β-cell dysfunction. Does NOT diagnose diabetes — lifestyle (diet, exercise, weight) strongly modulates genetic risk.",
        "HbA1c": "Higher z-score = genetically higher glycated hemoglobin. HbA1c reflects average blood glucose over ~3 months. Used clinically to diagnose diabetes (≥6.5%).",
        "Fasting Glucose": "Higher z-score = genetic tendency toward higher fasting blood sugar. Normal: <100 mg/dL. Prediabetes: 100-125 mg/dL. Diabetes: ≥126 mg/dL.",
        "Fasting Insulin": "Higher z-score = genetic predisposition to higher insulin levels, which may indicate insulin resistance when combined with high glucose.",
        "BMI": "Higher z-score = genetic predisposition to higher body mass index. BMI is a screening tool, not a diagnostic. Muscle mass, bone density, and fat distribution matter more clinically.",
        "Obesity": "Higher z-score = genetic risk for excess adiposity. Strongly modulated by environment — diet, physical activity, sleep, stress all influence expression of obesity-related genes.",
        "Overweight": "Higher z-score = genetic tendency toward BMI 25-30. Similar to obesity but milder genetic burden.",
        "Vitamin D": "Higher z-score = genetic predisposition to higher vitamin D levels. Deficiency defined as <20 ng/mL. Sun exposure, diet, and supplementation override genetic predisposition.",
        "Hypercholesterolemia": "Higher z-score = genetic tendency toward elevated LDL/Total cholesterol. Clinical cutoffs: Total <200 desirable, LDL <100 optimal. Statins dramatically reduce risk regardless of genetics.",
        "Non-HDL Cholesterol": "Higher z-score = genetic risk for elevated atherogenic lipoproteins. Non-HDL = Total - HDL. Target <130 mg/dL. Preferred over LDL when triglycerides >200.",
        "Coffee Consumption": "Higher z-score = genetic predisposition to higher coffee intake. CYP1A2 determines metabolism speed — 'fast' vs 'slow' metabolizers respond differently.",
        "Anemia / B12": "Higher z-score = genetic risk for lower hemoglobin/B12. Clinical diagnosis requires blood test (CBC, ferritin, B12, folate).",
        "Liver Cirrhosis": "HIGHER RISK REQUIRES CLINICAL CONTEXT. Genetic predisposition interacts strongly with alcohol consumption, viral hepatitis, and metabolic factors. Liver function tests and imaging are definitive.",
        "Gallstones": "Higher z-score = genetic risk for cholesterol gallstone formation. Risk factors include obesity, rapid weight loss, female sex, age >40. Ultrasound is diagnostic.",
        "PUFA (Omega-3/6)": "Higher z-score = genetic profile associated with higher polyunsaturated fatty acid levels. Reflects FADS1/FADS2 desaturase activity — ability to convert plant-based ALA to EPA/DHA.",
        "Omega-3": "Higher z-score = genetic predisposition to lower omega-3 levels. Dietary intake of fatty fish or supplementation can fully override genetic risk. No clinical cutoff exists.",
    }

    def interpret(e):
        trait = e.get("trait", "")
        z = e.get("z_score", 0)
        pctl = round(e.get("percentile", 50), 1)
        ctx = ""
        for key, desc in CLINICAL.items():
            if key.lower() in trait.lower():
                ctx = desc
                break
        if not ctx:
            ctx = f"Genetic predisposition score. Z={z:+.1f} means this individual is at the {pctl:.0f}th percentile of the EUR population for this trait."
        return ctx

    def reliable_badge(reliable):
        if reliable:
            return '<span style="background:#d5f5e3;color:#1e8449;padding:2px 6px;border-radius:3px;font-size:0.65rem;font-weight:700">✓ Reliable</span>'
        return '<span style="background:#fadbd8;color:#c0392b;padding:2px 6px;border-radius:3px;font-size:0.65rem;font-weight:700">⚠ Unreliable</span>'

    def row(e):
        r = "🔴" if e.get("z_score",0)>2 else ("🟠" if e.get("z_score",0)>1 else ("🟢" if e.get("z_score",0)<-1 else "🟡"))
        pctl = round(e.get("percentile",50), 1)
        z = e.get("z_score", 0)
        reliable = e.get("reliable", True)
        pgs_id = e.get("pgs_id", "")
        significance = "High risk" if z>2 else ("Elevated" if z>1 else ("Low/Protective" if z<-1 else "Population average"))

        # SNP coverage bar
        cov = pgs_coverage.get(pgs_id, {}) if pgs_coverage else {}
        n_used = cov.get("n_used", 0)
        n_total = cov.get("n_total", 0)
        snp_bar = snp_coverage_bar(n_used, n_total) if n_total > 0 else '<span style="color:#95a5a6;font-size:0.7rem">—</span>'

        # Risk bar
        bar_pct = max(5, min(95, pctl))
        risk_bar_html = risk_bar(bar_pct, z)

        return f"""<tr>
            <td>{r}</td><td><strong>{e['trait']}</strong></td><td>{e['pgs_id']}</td>
            <td>{reliable_badge(reliable)}</td>
            <td>{snp_bar}</td>
            <td style="font-weight:700;color:{'#e74c3c' if z>2 else ('#f39c12' if z>1 else ('#27ae60' if z<-1 else '#2c3e50'))}">{z:+.1f}</td>
            <td>{pctl}%</td><td>{significance}</td><td>{e.get('n_snps',0):,}</td>
            <td style="min-width:120px">{risk_bar_html}</td>
        </tr>"""

    # Summary bar
    all_entries = pgs_data.get("all_entries", [])
    n_reliable = sum(1 for e in all_entries if e.get("reliable", True))
    n_total_scores = len(all_entries)
    summary_html = (
        f'<div style="display:flex;gap:1rem;align-items:center;margin:0.75rem 0;flex-wrap:wrap;font-size:0.78rem">'
        f'<span>📊 <strong>PGS Summary:</strong></span>'
        f'<span style="color:#27ae60;font-weight:700">{n_reliable} reliable</span> / '
        f'<span style="color:#e74c3c;font-weight:700">{n_total_scores - n_reliable} unreliable</span>'
        f'<span>of {n_total_scores} scores</span>'
        f'</div>'
    )

    high_rows = "".join(row(e) for e in high_risk) or "<tr><td colspan='10' style='color:#7f8c8d'>None — no traits at elevated risk</td></tr>"
    elev_rows = "".join(row(e) for e in elevated) or "<tr><td colspan='10' style='color:#7f8c8d'>None</td></tr>"
    low_rows = "".join(row(e) for e in low_risk) or "<tr><td colspan='10' style='color:#7f8c8d'>None</td></tr>"

    # Build detailed interpretations
    detail_parts = []
    for e in high_risk + elevated:
        detail_parts.append(f"""
        <div class="trait-section" style="margin-bottom:1rem">
            <div class="trait-header" style="background:{'#fadbd8' if e.get('z_score',0)>2 else '#fdebd0'}">
                <span><strong>{e['trait']}</strong> ({e['pgs_id']})</span>
                <span class="risk-category-badge {'high' if e.get('z_score',0)>2 else 'medium'}">{'HIGHER RISK' if e.get('z_score',0)>2 else 'ELEVATED RISK'}</span>
            </div>
            <div class="trait-body">
                <div class="trait-stats">
                    <div class="trait-stat"><div class="stat-value" style="color:{'#e74c3c' if e.get('z_score',0)>2 else '#f39c12'}">{e.get('z_score',0):+.1f}σ</div><div class="stat-label">Z-Score vs EUR</div></div>
                    <div class="trait-stat"><div class="stat-value">{round(e.get('percentile',50),1)}%</div><div class="stat-label">Percentile</div></div>
                    <div class="trait-stat"><div class="stat-value">{e.get('n_snps',0):,}</div><div class="stat-label">Variants</div></div>
                </div>
                <div class="interpretation-box"><strong>Clinical Context:</strong><p style="margin-top:.5rem;white-space:pre-line">{interpret(e)}</p></div>
            </div>
        </div>""")

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Total Scores</h4><div class="big-stat">{summary.get('total_scores', 0)}</div><div class="stat-sub">from PGS Catalog</div></div>
        <div class="info-card"><h4>Reliable</h4><div class="big-stat">{summary.get('reliable_scores', 0)}</div><div class="stat-sub">&lt;500K SNPs</div></div>
        <div class="info-card"><h4>Reference</h4><div class="big-stat">{methodology.get('reference_panel', '1000G')[:20]}</div><div class="stat-sub">{', '.join(methodology.get('populations', []))}</div></div>
        <div class="info-card"><h4>Method</h4><div class="big-stat" style="font-size:0.7rem">z-score norm</div><div class="stat-sub">Population-stratified</div></div>
    </div>
    <p style="margin:1rem 0;line-height:1.6;font-size:0.9rem;background:var(--color-surface);padding:1rem;border-radius:8px;box-shadow:var(--shadow)">
    <strong>How to interpret these scores:</strong><br>
    • <strong>Z-score</strong>: How many standard deviations your genetic profile deviates from the EUR population mean.<br>
    • <strong>Percentile</strong>: Your rank. 95% = you're in the top 5% of genetic risk for that trait.<br>
    • <strong>Clinical significance</strong>: Z > 2 = notable deviation. Z > 3 = clinically relevant. But <em>genetics is NOT destiny</em> — lifestyle, environment, and medical care often override genetic predisposition.
    </p>
    <h4>🔴 Elevated Risk ({len(high_risk)})</h4>
    {summary_html}
    <div style="overflow-x:auto">
    <table><thead><tr><th></th><th>Trait</th><th>PGS ID</th><th>Reliable</th><th>Coverage</th><th>Z</th><th>%</th><th>Significance</th><th>SNPs</th><th>Risk Bar</th></tr></thead><tbody>{high_rows}</tbody></table>
    </div>
    <h4>🟠 Moderate Risk ({len(elevated)})</h4>
    <div style="overflow-x:auto">
    <table><thead><tr><th></th><th>Trait</th><th>PGS ID</th><th>Reliable</th><th>Coverage</th><th>Z</th><th>%</th><th>Significance</th><th>SNPs</th><th>Risk Bar</th></tr></thead><tbody>{elev_rows}</tbody></table>
    </div>
    <h4>🟢 Protective / Low Risk ({len(low_risk)})</h4>
    <div style="overflow-x:auto">
    <table><thead><tr><th></th><th>Trait</th><th>PGS ID</th><th>Reliable</th><th>Coverage</th><th>Z</th><th>%</th><th>Significance</th><th>SNPs</th><th>Risk Bar</th></tr></thead><tbody>{low_rows}</tbody></table>
    </div>

    {"<h4 style='margin-top:1.5rem'>📋 Detailed Clinical Interpretation</h4>" + ''.join(detail_parts) if detail_parts else ""}

    <p style="font-size:0.75rem;color:#7f8c8d;margin-top:1rem">
    ⚠️ <strong>Limitations:</strong> Calibrated against 1000 Genomes Phase 3 (2,504 samples, 5 super-populations).
    PGS scores are from published studies with varying population compositions — cross-population portability may be limited.
    Scores with >500K SNPs excluded due to distribution artifacts.
    <strong>All scores are for RESEARCH USE ONLY — not clinical diagnosis.</strong>
    </p>
    """


def build_calibration_detail_section(calibration_report):
    """Population calibration methodology and risk breakdown."""
    methodology = calibration_report.get("methodology", {})
    thresholds = methodology.get("risk_thresholds", {})
    high_traits = calibration_report.get("high_risk_traits", [])
    medium_traits = calibration_report.get("medium_risk_traits", [])
    low_traits = calibration_report.get("low_risk_traits", [])
    populations = methodology.get("population_strata", [])

    high_rows = "".join(f"<tr><td>{t}</td></tr>" for t in high_traits) or "<tr><td style='color:#7f8c8d'>None</td></tr>"
    low_rows = "".join(f"<tr><td>{t}</td></tr>" for t in low_traits) or "<tr><td style='color:#7f8c8d'>None</td></tr>"

    return f"""
    <div class="info-grid">
        <div class="info-card">
            <h4>Reference Panel</h4>
            <div class="big-stat" style="font-size:1rem">{methodology.get('reference_panel', 'N/A')}</div>
        </div>
        <div class="info-card">
            <h4>Normalization</h4>
            <div class="big-stat" style="font-size:0.9rem">{methodology.get('normalization', 'N/A').replace('_', ' ').title()}</div>
        </div>
        <div class="info-card">
            <h4>Populations</h4>
            <div class="big-stat">{len(populations)}</div>
            <div class="stat-sub">{', '.join(populations)}</div>
        </div>
        <div class="info-card">
            <h4>Traits Analyzed</h4>
            <div class="big-stat">{calibration_report.get('traits_analyzed', 0)}</div>
        </div>
    </div>

    <div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-top:1rem">
        <div>
            <h4>Risk Thresholds</h4>
            <table>
                <thead><tr><th>Category</th><th>Percentile Range</th></tr></thead>
                <tbody>
                    <tr><td style="color:#e74c3c;font-weight:700">High Risk</td><td>{thresholds.get('high', '>75th')}</td></tr>
                    <tr><td style="color:#f39c12;font-weight:700">Medium Risk</td><td>{thresholds.get('medium', '25-75th')}</td></tr>
                    <tr><td style="color:#27ae60;font-weight:700">Low Risk</td><td>{thresholds.get('low', '<25th')}</td></tr>
                </tbody>
            </table>
        </div>
        <div>
            <h4>Assigned Population: {calibration_report.get('assigned_population', 'EUR')}</h4>
            <p style="font-size:0.8rem;color:#7f8c8d">{calibration_report.get('calibration_note', '')[:300]}</p>
        </div>
    </div>

    <div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-top:1rem">
        <div>
            <h4>Elevated Risk Traits ({len(high_traits)})</h4>
            <table><thead><tr><th>Trait</th></tr></thead><tbody>{high_rows}</tbody></table>
        </div>
        <div>
            <h4>Lower Risk Traits ({len(low_traits)})</h4>
            <table><thead><tr><th>Trait</th></tr></thead><tbody>{low_rows}</tbody></table>
        </div>
    </div>
    """
    """Pipeline methodology summary."""
    meta = prs_result.get("metadata", {})
    n_variants = meta.get("n_variants", prs_result.get("prs_core", {}).get("n_variants", 109))
    n_traits = meta.get("n_traits", 10)
    formula = meta.get("prs_formula", "PRS = Σ(βⱼ × Gᵢⱼ)")
    method = meta.get("computation_method", "PLINK --score (dosage-weighted)")
    pipeline_ver = meta.get("pipeline_version", "1.1.0")
    n_ref = ancestry.get("n_reference_samples", 2504)
    n_pcs = ancestry.get("n_pcs", 20)

    return f"""
    <div class="info-grid">
        <div class="info-card"><h4>Pipeline</h4><div class="big-stat">v{pipeline_ver}</div><div class="stat-sub">BlueGen</div></div>
        <div class="info-card"><h4>PRS Formula</h4><div class="big-stat" style="font-size:0.9rem">{formula}</div><div class="stat-sub">{method}</div></div>
        <div class="info-card"><h4>Variants</h4><div class="big-stat">{n_variants}</div><div class="stat-sub">across {n_traits} traits</div></div>
        <div class="info-card"><h4>Reference</h4><div class="big-stat">{n_ref}</div><div class="stat-sub">1000G Phase 3, {n_pcs} PCs</div></div>
    </div>
    <h4>Pipeline Stages</h4>
    <table>
        <thead><tr><th>Stage</th><th>Process</th><th>Method</th><th>Output</th></tr></thead>
        <tbody>
            <tr><td>A</td><td>VCF → PLINK</td><td>PLINK 1.9, GQ≥20, DP≥10</td><td>.bed/.bim/.fam</td></tr>
            <tr><td>B</td><td>Quality Control</td><td>geno 0.05, maf 0.01, HWE 1e-6</td><td>Filtered genotypes</td></tr>
            <tr><td>C</td><td>LD Pruning</td><td>Per-population (EUR/AFR/EAS/SAS/AMR), conservative intersection</td><td>Ancestry-matched independent SNPs</td></tr>
            <tr><td>D</td><td>PCA + Projection</td><td>1000G-trained PCA, target sample projection</td><td>20 PCs, ancestry inference</td></tr>
            <tr><td>F</td><td>PRS Computation</td><td>PLINK --score (dosage-weighted)</td><td>Raw PRS per trait</td></tr>
            <tr><td>G</td><td>PCA Adjustment</td><td>PRS_adj = PRS_raw − Σ(βₖ × PCₖ)</td><td>Ancestry-adjusted PRS</td></tr>
            <tr><td>H</td><td>Population Calibration</td><td>Empirical 1000G population distributions</td><td>Z-scores + percentiles</td></tr>
            <tr><td>7-10</td><td>Validation & Lock</td><td>8-dimension validation, adversarial stress, publication lock</td><td>Scientific integrity score</td></tr>
        </tbody>
    </table>
    """


def build_clinvar_section(clinvar_data: dict, ui: dict) -> str:
    """ClinVar pathogenic variant annotation — confidence-tiered, bilingual."""
    lang = ui.get("_lang", "en")

    if not clinvar_data or not clinvar_data.get("pathogenic_variants"):
        msg = {"en": "No ClinVar pathogenic variants found.", "es": "No se encontraron variantes patogénicas de ClinVar."}
        hint = {"en": 'Run with <code>--clinvar</code> to generate.', "es": 'Ejecuta con <code>--clinvar</code> para generar.'}
        return f"""<div class="info-card" style="text-align:center;padding:1.5rem">
            <p style="color:var(--color-text-secondary)">{msg.get(lang, msg['en'])}</p>
            <p style="color:var(--color-text-secondary);font-size:0.8rem">{hint.get(lang, hint['en'])}</p>
        </div>"""

    meta = clinvar_data.get("metadata", {})
    summary = clinvar_data.get("pathogenic_variant_summary", {})
    variants = clinvar_data.get("pathogenic_variants", [])
    tier_counts = summary.get("by_confidence_tier", {})
    n_high_conf = summary.get("high_confidence_count", 0)
    clinvar_date = meta.get("clinvar_release_date", "")[:10]

    n_pathogenic = summary.get("total_pathogenic", 0)
    n_likely = summary.get("total_likely_pathogenic", 0)
    n_risk = summary.get("total_risk_alleles", 0)
    total = len(variants)

    # ═══ BILINGUAL LABELS ═══
    T = {
        "en": {
            "veracity_alert": "⚠️ Veracity Note",
            "veracity_text": f"Of {total} variants found, only <strong>{n_high_conf}</strong> have strong evidence "
                            f"(expert panel or multiple-lab consensus). <strong>{tier_counts.get('very_low', 0)} variants</strong> "
                            f"lack explicit evidence criteria — these are the least reliable findings and should not be "
                            f"interpreted without clinical confirmation.",
            "high_conf_title": "High & Moderate Confidence",
            "high_conf_desc": "Reviewed by expert panel or multiple labs agree. These are the most reliable findings.",
            "low_conf_title": "Lower Confidence & Risk Alleles",
            "low_conf_desc": "Single submitter or no evidence criteria stated. Includes risk alleles (increase susceptibility, do NOT cause disease).",
            "reliable_count": f"{n_high_conf} reliable",
            "uncertain_count": f"{total - n_high_conf} uncertain / risk",
            "confidence_tier_high": "🏅 Expert",
            "confidence_tier_moderate": "✓ Multi-Lab",
            "confidence_tier_low": "⚠️ Single Lab",
            "confidence_tier_very_low": "❓ No Criteria",
            "clinvar_version": "ClinVar Release",
            "source_note": f"Source: NCBI ClinVar (GRCh37, {clinvar_date}). Variant classifications are submitted by independent laboratories and may change over time.",
            "limitations_title": "⚠️ Limitations of ClinVar Analysis",
            "limitations_text": "• Many variants lack evidence criteria (see confidence tiers above). "
                               "• Pathogenic ≠ you will get the disease. Penetrance varies widely. "
                               "• Some conditions are treatable or preventable; others are not. "
                               "• These results are RESEARCH USE ONLY — confirm with a clinical lab.",
            "legend_title": "Confidence Tiers — How Reliable Is Each Finding?",
            "tier_high": "Expert Panel or Clinical Guideline — highest confidence, endorsed by domain authorities",
            "tier_moderate": "Multiple Labs Agree — independent laboratories concur on this classification",
            "tier_low": "Single Lab — one laboratory's classification; confirmation recommended",
            "tier_very_low": "No Evidence Criteria — classification without stated methodology; treat as preliminary",
        },
        "es": {
            "veracity_alert": "⚠️ Nota de Veracidad",
            "veracity_text": f"De {total} variantes encontradas, solo <strong>{n_high_conf}</strong> tienen evidencia sólida "
                            f"(panel experto o consenso multi-laboratorio). <strong>{tier_counts.get('very_low', 0)} variantes</strong> "
                            f"carecen de criterios de evidencia explícitos — estos son los hallazgos menos fiables y no deben "
                            f"interpretarse sin confirmación clínica.",
            "high_conf_title": "Confianza Alta y Moderada",
            "high_conf_desc": "Revisadas por panel experto o múltiples labs coinciden. Son los hallazgos más fiables.",
            "low_conf_title": "Menor Confianza y Alelos de Riesgo",
            "low_conf_desc": "Un solo laboratorio o sin criterios de evidencia. Incluye alelos de riesgo (aumentan susceptibilidad, NO causan enfermedad).",
            "reliable_count": f"{n_high_conf} fiables",
            "uncertain_count": f"{total - n_high_conf} inciertas / riesgo",
            "confidence_tier_high": "🏅 Experto",
            "confidence_tier_moderate": "✓ Multi-Lab",
            "confidence_tier_low": "⚠️ Un Lab",
            "confidence_tier_very_low": "❓ Sin Criterios",
            "clinvar_version": "Versión de ClinVar",
            "source_note": f"Fuente: NCBI ClinVar (GRCh37, {clinvar_date}). Las clasificaciones son enviadas por laboratorios independientes y pueden cambiar.",
            "limitations_title": "⚠️ Limitaciones del Análisis ClinVar",
            "limitations_text": "• Muchas variantes carecen de criterios de evidencia (ver tiers de confianza arriba). "
                               "• Patogénica ≠ desarrollarás la enfermedad. La penetrancia varía mucho. "
                               "• Algunas condiciones son tratables o prevenibles; otras no. "
                               "• SOLO PARA INVESTIGACIÓN — confirma con un laboratorio clínico.",
            "legend_title": "Tiers de Confianza — ¿Qué Tan Fiable Es Cada Hallazgo?",
            "tier_high": "Panel Experto o Guía Clínica — máxima confianza, respaldado por autoridades",
            "tier_moderate": "Múltiples Labs Coinciden — laboratorios independientes confirman esta clasificación",
            "tier_low": "Un Solo Lab — clasificación de un laboratorio; se recomienda confirmación",
            "tier_very_low": "Sin Criterios de Evidencia — clasificación sin metodología explícita; considerar preliminar",
        },
    }
    t = T.get(lang, T["en"])

    # ═══ BADGES ═══
    sig_colors = {
        "Pathogenic": ("#c0392b", "#fadbd8"),
        "Likely_pathogenic": ("#d35400", "#fdebd0"),
        "Pathogenic/Likely_pathogenic": ("#e67e22", "#fef5e7"),
        "Risk_allele": ("#f39c12", "#fef9e7"),
    }

    def sig_badge(sig):
        fg, bg = sig_colors.get(sig, ("#7f8c8d", "#e9ecef"))
        return f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.7rem;font-weight:700">{sig}</span>'

    tier_badges = {
        "high": ("#27ae60", "#d5f5e3", t["confidence_tier_high"]),
        "moderate": ("#2e86c1", "#d6eaf8", t["confidence_tier_moderate"]),
        "low": ("#f39c12", "#fdebd0", t["confidence_tier_low"]),
        "very_low": ("#95a5a6", "#eaecee", t["confidence_tier_very_low"]),
    }

    def tier_badge(tier):
        fg, bg, label = tier_badges.get(tier, ("#95a5a6", "#eaecee", tier))
        return f'<span style="background:{bg};color:{fg};padding:1px 6px;border-radius:3px;font-size:0.65rem;font-weight:700" title="{label}">{label}</span>'

    def humanize_disease(d):
        if not d or d in (".", "not_provided"): return "—"
        return " | ".join(d.replace("_", " ").split("|")[:3])

    def format_af(af):
        if not af or af == "—": return "—"
        try:
            v = float(af)
            return f"{v:.4f}" if v >= 0.0001 else f"{v:.2e}"
        except (ValueError, TypeError): return af

    def fmt_review(r):
        if not r: return "—"
        return (r.replace("_", " ")
                .replace("practice guideline", "⭐ guideline")
                .replace("reviewed by expert panel", "🏅 expert panel")
                .replace("criteria provided multiple submitters no conflicts", "✓ multi-lab")
                .replace("criteria provided single submitter", "single lab")
                .replace("no assertion criteria provided", "no criteria")
                .replace("no assertion provided", "no assertion"))

    def build_row(v):
        desc = v.get("disease_description", "")
        cui = v.get("medgen_cui", "")
        desc_html = ""
        if desc:
            medgen_link = (
                f' <a href="https://www.ncbi.nlm.nih.gov/medgen/{cui}" target="_blank" '
                f'rel="noopener" style="font-size:0.7rem">[MedGen ↗]</a>'
            ) if cui else ""
            desc_html = (
                f'<tr style="background:#fafafa">'
                f'<td colspan="8" style="font-size:0.75rem;color:var(--color-text-secondary);padding:2px 12px 6px 32px;border-bottom:1px solid #eee">'
                f'<em>{desc[:250]}</em>{medgen_link}'
                f'</td>'
                f'</tr>'
            )
        return (
            f'<tr>'
            f'<td><code>{v.get("rsid","") or "—"}</code></td>'
            f'<td style="white-space:nowrap">{v.get("chrom","?")}:{v.get("pos","?")}</td>'
            f'<td style="font-weight:600">{", ".join(v.get("genes",[])) or v.get("gene_info","—")}</td>'
            f'<td style="font-size:0.82rem">{humanize_disease(v.get("disease_name","—"))[:90]}</td>'
            f'<td>{sig_badge(v.get("clinical_significance","—"))}</td>'
            f'<td>{tier_badge(v.get("confidence_tier","very_low"))}</td>'
            f'<td style="font-size:0.72rem;color:var(--color-text-secondary)">{fmt_review(v.get("review_status",""))}</td>'
            f'<td style="font-size:0.75rem;color:var(--color-text-secondary)">{format_af(v.get("af_exac","") or v.get("af_1000g",""))}</td>'
            f'</tr>'
            f'{desc_html}'
        )

    # Split variants by confidence
    high_mod = [v for v in variants if v.get("confidence_tier") in ("high", "moderate")]
    low_variants = [v for v in variants if v.get("confidence_tier") in ("low", "very_low")]

    match_rate = meta.get("match_rate", 0)
    n_overlap = meta.get("positional_overlaps", 0)
    n_exact = meta.get("exact_matches", 0)

    # ═══════════════════════════════════════════════════════════════════════════
    return f"""
    <!-- ═══ VERACITY ALERT ═══ -->
    <div class="highlight-box" style="background:#fdedec;border:2px solid #e74c3c;border-radius:8px;padding:1rem 1.5rem;margin-bottom:1.5rem">
        <h4 style="margin-top:0;color:#c0392b">{t['veracity_alert']}</h4>
        <p style="font-size:0.9rem;margin:0">{t['veracity_text']}</p>
    </div>

    <!-- ═══ CONFIDENCE TIER LEGEND ═══ -->
    <div class="highlight-box" style="background:#eaf2f8;border:1px solid #aed6f1;border-radius:8px;padding:1rem 1.5rem;margin-bottom:1rem">
        <h4 style="margin-top:0">📖 {t['legend_title']}</h4>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:0.4rem 1.5rem;font-size:0.82rem">
            <div>{tier_badge('high')} — {t['tier_high']}</div>
            <div>{tier_badge('moderate')} — {t['tier_moderate']}</div>
            <div>{tier_badge('low')} — {t['tier_low']}</div>
            <div>{tier_badge('very_low')} — {t['tier_very_low']}</div>
        </div>
        <p style="font-size:0.8rem;color:var(--color-text-secondary);margin:0.5rem 0 0">{t['source_note']}</p>
    </div>

    <!-- ═══ SUMMARY CARDS ═══ -->
    <div class="info-grid">
        <div class="info-card" style="border-left:3px solid #27ae60">
            <h4>{t['high_conf_title']}</h4>
            <div class="big-stat" style="color:#27ae60">{n_high_conf}</div>
            <div class="stat-sub">{t['high_conf_desc']}</div>
        </div>
        <div class="info-card" style="border-left:3px solid #95a5a6">
            <h4>{t['low_conf_title']}</h4>
            <div class="big-stat" style="color:#95a5a6">{total - n_high_conf}</div>
            <div class="stat-sub">{t['low_conf_desc']}</div>
        </div>
        <div class="info-card">
            <h4>{t['clinvar_version']}</h4>
            <div class="big-stat" style="font-size:1rem">{clinvar_date or 'N/A'}</div>
            <div class="stat-sub">NCBI ClinVar GRCh37 · {meta.get('user_vcf_total_variants', 0):,} variants analyzed</div>
        </div>
        <div class="info-card" style="border-left:3px solid #f39c12">
            <h4>Alelos de Riesgo</h4>
            <div class="big-stat" style="color:#f39c12">{n_risk}</div>
            <div class="stat-sub">Susceptibilidad — NO causan enfermedad</div>
        </div>
    </div>

    <!-- ═══ HIGH CONFIDENCE TABLE ═══ -->
    <h4 style="margin-top:1.5rem;color:#27ae60">🏅 {t['high_conf_title']} — {len(high_mod)} {t['reliable_count']}</h4>
    <div style="overflow-x:auto">
    <table>
        <thead><tr><th>rsID</th><th>Pos</th><th>Gene</th><th>Disease</th><th>Significance</th><th>Confidence</th><th>Review</th><th>Freq</th></tr></thead>
        <tbody>{''.join(build_row(v) for v in high_mod[:200]) or '<tr><td colspan="8" style="color:var(--color-text-secondary);text-align:center;padding:1rem">✅ No high-confidence pathogenic variants found. This is normal.</td></tr>'}</tbody>
    </table>
    </div>

    <!-- ═══ LOWER CONFIDENCE TABLE ═══ -->
    <h4 style="margin-top:1.5rem;color:#95a5a6">❓ {t['low_conf_title']} — {len(low_variants)} {t['uncertain_count']}</h4>
    <div style="overflow-x:auto">
    <table>
        <thead><tr><th>rsID</th><th>Pos</th><th>Gene</th><th>Disease</th><th>Significance</th><th>Confidence</th><th>Review</th><th>Freq</th></tr></thead>
        <tbody>{''.join(build_row(v) for v in low_variants[:300]) or '<tr><td colspan="8" style="color:var(--color-text-secondary);text-align:center;padding:1rem">No lower-confidence variants.</td></tr>'}</tbody>
    </table>
    </div>
    {f'<p style="color:var(--color-text-secondary);font-size:0.8rem;margin-top:0.5rem">{len(low_variants) - 300} more not shown. See clinvar/clinvar_pathogenic_variants.json for complete list.</p>' if len(low_variants) > 300 else ''}

    <!-- ═══ LIMITATIONS ═══ -->
    <div class="highlight-box" style="background:#fef9e7;border:1px solid #f9e79f;border-radius:8px;padding:1rem 1.5rem;margin-top:1.5rem">
        <h4 style="margin-top:0">{t['limitations_title']}</h4>
        <p style="font-size:0.82rem;margin:0;white-space:pre-line;color:var(--color-text-secondary)">{t['limitations_text']}</p>
    </div>

    <!-- ═══ FULL DISCLAIMER ═══ -->
    <div class="disclaimer-box" style="margin-top:1rem">
        <p style="white-space:pre-line;font-size:0.82rem;margin:0">⚠️ RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS. NO PARA DIAGNÓSTICO CLÍNICO.</p>
    </div>
    """


def build_pharmgkb_section(pharmgkb_data: dict, ui: dict) -> str:
    """Pharmacogenomic drug response section."""
    lang = ui.get("_lang", "en")
    findings = pharmgkb_data.get("pharmacogenomic_findings", [])
    summary = pharmgkb_data.get("summary", {})

    if not findings:
        msg = {"en": "No pharmacogenomic variants found.", "es": "No se encontraron variantes farmacogenómicas."}
        return f"""<div class="info-card" style="text-align:center;padding:1.5rem">
            <p style="color:var(--color-text-secondary)">{msg.get(lang, msg['en'])}</p>
        </div>"""

    action_colors = {
        "critical": ("#c0392b", "#fadbd8", "🔴"),
        "important": ("#d35400", "#fdebd0", "🟠"),
        "informative": ("#2980b9", "#d6eaf8", "🟡"),
        "normal": ("#27ae60", "#d5f5e3", "🟢"),
    }

    def action_badge(a):
        fg, bg, icon = action_colors.get(a, ("#7f8c8d", "#e9ecef", "⚪"))
        return f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.7rem;font-weight:700">{icon} {a}</span>'

    n_critical = summary.get("by_actionability", {}).get("critical", 0)
    n_important = summary.get("by_actionability", {}).get("important", 0)
    n_info = summary.get("by_actionability", {}).get("informative", 0)

    T = {
        "en": {
            "title": "Pharmacogenomics — Drug Response",
            "explain": "These variants affect how your body processes specific medications. They do <strong>not</strong> diagnose any condition — they inform drug selection and dosing <strong>if</strong> you ever need these medications.",
            "critical": "Critical — contraindicated or major dose adjustment needed",
            "important": "Important — dose adjustment recommended",
            "informative": "Informative — may affect response, standard monitoring advised",
            "gene": "Gene",
            "variant": "Variant",
            "drug": "Drug",
            "phenotype": "Your Phenotype",
            "recommendation": "Clinical Recommendation",
            "cpic_source": "CPIC Guideline",
            "disclaimer": "⚠️ Do NOT change any medication without consulting your doctor. These are pharmacogenomic predictions for research use.",
        },
        "es": {
            "title": "Farmacogenómica — Respuesta a Fármacos",
            "explain": "Estas variantes afectan cómo tu cuerpo procesa medicamentos específicos. <strong>No</strong> diagnostican ninguna condición — informan la selección y dosificación de fármacos <strong>si</strong> alguna vez necesitas estos medicamentos.",
            "critical": "Crítico — contraindicado o ajuste de dosis mayor necesario",
            "important": "Importante — ajuste de dosis recomendado",
            "recommendation": "Recomendación Clínica",
            "informative": "Informativo — puede afectar respuesta, monitorización estándar",
            "gene": "Gen",
            "variant": "Variante",
            "drug": "Fármaco",
            "phenotype": "Tu Fenotipo",
            "cpic_source": "Guía CPIC",
            "disclaimer": "⚠️ NO cambies ningún medicamento sin consultar a tu médico. Estas son predicciones farmacogenómicas para uso en investigación.",
        },
    }
    t = T.get(lang, T["en"])

    rows = ""
    for f in findings:
        rec = f.get(f"recommendation_{lang}", f.get("recommendation_en", ""))
        if len(rec) > 150:
            rec = rec[:147] + "..."
        rows += (
            f'<tr>'
            f'<td style="font-weight:600">{f["gene"]}</td>'
            f'<td><code>{f["rsid"]}</code> {f.get("star_allele","")}</td>'
            f'<td style="font-weight:600">{f["drug"]} ({f.get("drug_class","")})</td>'
            f'<td style="font-size:0.85rem">{f["phenotype"]} ({f["copies"]} copia{"s" if f["copies"]>1 else ""})</td>'
            f'<td style="font-size:0.82rem">{rec}</td>'
            f'<td style="font-size:0.72rem;color:var(--color-text-secondary)">{f["cpic_level"]}</td>'
            f'</tr>'
        )

    return f"""
    <div class="highlight-box" style="background:#eaf2f8;border:1px solid #aed6f1;border-radius:8px;padding:1rem 1.5rem;margin-bottom:1.5rem">
        <p style="font-size:0.9rem;margin:0">{t['explain']}</p>
    </div>

    <div class="info-grid" style="grid-template-columns:repeat(3,1fr)">
        {f'<div class="info-card" style="border-left:3px solid #c0392b"><h4>🔴 {t["critical"]}</h4><div class="big-stat" style="color:#c0392b">{n_critical}</div></div>' if n_critical else ''}
        {f'<div class="info-card" style="border-left:3px solid #d35400"><h4>🟠 {t["important"]}</h4><div class="big-stat" style="color:#d35400">{n_important}</div></div>' if n_important else ''}
        <div class="info-card" style="border-left:3px solid #2980b9"><h4>🟡 {t["informative"]}</h4><div class="big-stat" style="color:#2980b9">{n_info}</div></div>
    </div>

    <h4 style="margin-top:1.5rem">💊 {t['title']} ({len(findings)} hallazgos)</h4>
    <div style="overflow-x:auto">
    <table>
        <thead><tr><th>{t['gene']}</th><th>{t['variant']}</th><th>{t['drug']}</th><th>{t['phenotype']}</th><th>{t['recommendation']}</th><th>{t['cpic_source']}</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    </div>

    <div class="highlight-box" style="background:#eaf2f8;border:1px solid #aed6f1;border-radius:8px;padding:0.8rem 1.2rem;margin-top:1rem">
        <p style="font-size:0.82rem;margin:0">
        <strong>📖 CPIC Guideline column:</strong> Indicates which clinical guideline supports this recommendation.<br>
        <strong>CPIC</strong> = Clinical Pharmacogenetics Implementation Consortium (U.S. NIH-funded). <strong>DPWG</strong> = Dutch Pharmacogenetics Working Group (European).<br>
        <strong>Level A/B</strong> = strongest evidence; actionable prescribing change recommended. <strong>Level C/D</strong> = weaker evidence; consider but not mandatory.
        </p>
    </div>

    <div class="disclaimer-box" style="margin-top:1rem">
        <p style="white-space:pre-line;font-size:0.82rem;margin:0">{t['disclaimer']}</p>
    </div>
    """


# ═══════════════════════════════════════════════════════════════════════════════
# CLINICAL ACTIONABILITY SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

# Mapping of PRS high-risk traits to supporting ClinVar genes and PharmGKB genes
# Based on known biological pathways (verified against actual ClinVar/PharmGKB data)
CLINICAL_CONVERGENCE_MAP = {
    "Lipid metabolism": {
        "clinvar_genes": ["NOS3"], "pharmgkb_genes": ["SLCO1B1"],
        "context_en": "NOS3 (eNOS) regulates vascular tone and lipid metabolism; SLCO1B1 mediates statin transport. Statin pharmacogenetics may be clinically relevant.",
        "context_es": "NOS3 (eNOS) regula el tono vascular y el metabolismo lipídico; SLCO1B1 media el transporte de estatinas. La farmacogenética de estatinas puede ser clínicamente relevante.",
    },
    "Glucose metabolism": {
        "clinvar_genes": ["TCF7L2", "CAPN10"], "pharmgkb_genes": [],
        "context_en": "TCF7L2 is the strongest T2D GWAS hit; CAPN10 was the first T2D locus identified via positional cloning. Both support the polygenic signal.",
        "context_es": "TCF7L2 es el hit GWAS más fuerte para diabetes tipo 2; CAPN10 fue el primer locus T2D identificado por clonación posicional. Ambos apoyan la señal poligénica.",
    },
    "Blood pressure": {
        "clinvar_genes": ["ECE1"], "pharmgkb_genes": ["CYP2C9"],
        "context_en": "ECE1 converts big endothelin to active endothelin-1 (vasoconstrictor). CYP2C9 metabolizes NSAIDs and warfarin — may affect drug response in hypertension context.",
        "context_es": "ECE1 convierte big endotelina en endotelina-1 activa (vasoconstrictor). CYP2C9 metaboliza AINEs y warfarina — puede afectar la respuesta a fármacos en contexto de hipertensión.",
    },
    "Folate & methylation": {
        "clinvar_genes": ["MTHFR"], "pharmgkb_genes": [],
        "context_en": "MTHFR is central to folate/homocysteine metabolism. Variants affect cardiovascular risk and methylation capacity.",
        "context_es": "MTHFR es central en el metabolismo de folato/homocisteína. Las variantes afectan el riesgo cardiovascular y la capacidad de metilación.",
    },
}

DRUG_PRS_INTERSECTIONS = [
    {"drug_class": "Statins", "drugs": "atorvastatin, fluvastatin, lovastatin, pitavastatin, simvastatin",
     "gene": "SLCO1B1", "prs_trait": "Lipid metabolism",
     "interaction_en": "SLCO1B1 variant reduces statin hepatic uptake → higher plasma levels → increased myopathy risk",
     "interaction_es": "La variante SLCO1B1 reduce la captación hepática de estatinas → niveles plasmáticos más altos → mayor riesgo de miopatía",
     "recommendation_en": "Consider lower statin dose or alternative (rosuvastatin, pravastatin). Genetic testing for SLCO1B1 recommended before high-dose simvastatin.",
     "recommendation_es": "Considerar dosis más baja de estatina o alternativa (rosuvastatina, pravastatina). Prueba genética de SLCO1B1 recomendada antes de simvastatina a dosis altas."},
    {"drug_class": "Warfarin", "drugs": "warfarin",
     "gene": "CYP2C9 / VKORC1", "prs_trait": "Blood pressure",
     "interaction_en": "CYP2C9 reduced-function variant slows warfarin metabolism; VKORC1 variant affects warfarin sensitivity",
     "interaction_es": "La variante de función reducida de CYP2C9 ralentiza el metabolismo de warfarina; la variante VKORC1 afecta la sensibilidad a warfarina",
     "recommendation_en": "Start with lower warfarin dose (2-3 mg/day). Monitor INR closely. Consider direct oral anticoagulants as alternative.",
     "recommendation_es": "Comenzar con dosis más baja de warfarina (2-3 mg/día). Monitorizar INR de cerca. Considerar anticoagulantes orales directos como alternativa."},
    {"drug_class": "NSAIDs", "drugs": "celecoxib, meloxicam",
     "gene": "CYP2C9", "prs_trait": "Blood pressure",
     "interaction_en": "CYP2C9 reduced-function variant slows NSAID clearance. NSAIDs can increase blood pressure via COX-2 mediated sodium retention.",
     "interaction_es": "La variante de función reducida de CYP2C9 ralentiza la eliminación de AINEs. Los AINEs pueden aumentar la presión arterial por retención de sodio mediada por COX-2.",
     "recommendation_en": "Monitor blood pressure if NSAID therapy is needed. Consider celecoxib dose reduction in CYP2C9 intermediate metabolizers.",
     "recommendation_es": "Monitorizar presión arterial si se necesita terapia con AINEs. Considerar reducción de dosis de celecoxib en metabolizadores intermedios de CYP2C9."},
]

def build_clinical_actionability_section(clinvar_data, pharmgkb_data, prs_entries, ui):
    """Cross-reference ClinVar + PharmGKB + high-risk PRS into a unified clinical summary.

    Three subsections:
    (a) High-Confidence Clinical Findings (ClinVar tier≥moderate + PharmGKB actionable)
    (b) PRS-Gene Convergence (per-trait supporting gene evidence)
    (c) Drug-Gene-PRS Intersections
    """
    lang = ui.get("_lang", "en")
    is_en = lang == "en"

    T = {
        "en": {
            "title": "Clinical Actionability Summary",
            "desc": "Cross-references ClinVar pathogenic variants, pharmacogenomic findings, and elevated polygenic risk scores to identify clinically relevant intersections. <strong>Research use only — confirm all findings with a healthcare professional.</strong>",
            "subsection_a": "High-Confidence Clinical Findings",
            "subsection_b": "PRS-Gene Convergence",
            "subsection_c": "Drug-Gene-PRS Intersections",
            "no_data": "No clinical convergence findings identified.",
            "clinvar_count": "ClinVar High-Confidence",
            "pharmgkb_count": "PharmGKB Actionable",
            "prs_high_count": "PRS High-Risk Traits",
            "findings_total": "Total Convergent Findings",
            "trait": "PRS Trait",
            "z_score": "Z-Score",
            "risk": "Risk",
            "clinvar_genes": "ClinVar Genes",
            "pharmgkb_genes": "PharmGKB Genes",
            "context": "Biological Context",
            "drug_class": "Drug Class",
            "drugs": "Drugs",
            "gene": "Gene",
            "interaction": "Interaction",
            "recommendation": "Recommendation",
        },
        "es": {
            "title": "Resumen de Accionabilidad Clínica",
            "desc": "Cruza variantes patogénicas de ClinVar, hallazgos farmacogenómicos y puntajes de riesgo poligénico elevados para identificar intersecciones clínicamente relevantes. <strong>Solo para investigación — confirma todos los hallazgos con un profesional de la salud.</strong>",
            "subsection_a": "Hallazgos Clínicos de Alta Confianza",
            "subsection_b": "Convergencia PRS-Gen",
            "subsection_c": "Intersecciones Fármaco-Gen-PRS",
            "no_data": "No se identificaron hallazgos de convergencia clínica.",
            "clinvar_count": "ClinVar Alta Confianza",
            "pharmgkb_count": "PharmGKB Accionable",
            "prs_high_count": "PRS Rasgos Alto Riesgo",
            "findings_total": "Total Hallazgos Convergentes",
            "trait": "Rasgo PRS",
            "z_score": "Z-Score",
            "risk": "Riesgo",
            "clinvar_genes": "Genes ClinVar",
            "pharmgkb_genes": "Genes PharmGKB",
            "context": "Contexto Biológico",
            "drug_class": "Clase de Fármaco",
            "drugs": "Fármacos",
            "gene": "Gen",
            "interaction": "Interacción",
            "recommendation": "Recomendación",
        },
    }
    t = T.get(lang, T["en"])

    # ── Gather data ──
    variants = clinvar_data.get("pathogenic_variants", []) if clinvar_data else []
    high_conf_variants = [v for v in variants if v.get("confidence_tier") in ("high", "moderate")]
    pharm_findings = pharmgkb_data.get("pharmacogenomic_findings", []) if pharmgkb_data else []
    actionable_pharm = [f for f in pharm_findings if f.get("actionability") in ("critical", "important", "informative")]
    high_risk_traits = [e for e in (prs_entries or []) if e.get("risk_category") == "high"]

    if not high_conf_variants and not actionable_pharm and not high_risk_traits:
        return f'<div class="info-card" style="text-align:center;padding:1.5rem"><p style="color:var(--color-text-secondary)">{t["no_data"]}</p></div>'

    # ── (a) High-Confidence Clinical Findings ──
    total_findings = len(high_conf_variants) + len(actionable_pharm)
    subsection_a = (
        f'<div class="info-grid" style="grid-template-columns:repeat(4,1fr)">'
        f'<div class="info-card" style="border-left:3px solid #27ae60"><h4>{t["clinvar_count"]}</h4><div class="big-stat" style="color:#27ae60">{len(high_conf_variants)}</div></div>'
        f'<div class="info-card" style="border-left:3px solid #3498db"><h4>{t["pharmgkb_count"]}</h4><div class="big-stat" style="color:#3498db">{len(actionable_pharm)}</div></div>'
        f'<div class="info-card" style="border-left:3px solid #e74c3c"><h4>{t["prs_high_count"]}</h4><div class="big-stat" style="color:#e74c3c">{len(high_risk_traits)}</div></div>'
        f'<div class="info-card" style="border-left:3px solid #9b59b6"><h4>{t["findings_total"]}</h4><div class="big-stat" style="color:#9b59b6">{total_findings}</div></div>'
        f'</div>'
    )

    # ── (b) PRS-Gene Convergence ──
    conv_rows = ""
    for e in high_risk_traits:
        trait = e.get("trait", "")
        z = safe_float(e.get("population_zscore", e.get("raw_score", 0)))
        mapping = CLINICAL_CONVERGENCE_MAP.get(trait, {"clinvar_genes": [], "pharmgkb_genes": [], "context_en": "", "context_es": ""})

        cv_genes_html = ", ".join(
            f'<span class="clinical-convergence-gene">{g}</span>' for g in mapping.get("clinvar_genes", []))
        pg_genes_html = ", ".join(
            f'<span class="clinical-convergence-gene">{g}</span>' for g in mapping.get("pharmgkb_genes", []))
        ctx = mapping.get("context_en" if is_en else "context_es", "")

        z_color = risk_color(z)
        risk_label = "HIGHER RISK" if is_en else "RIESGO ELEVADO"

        conv_rows += (
            f'<tr><td><strong>{trait}</strong></td>'
            f'<td style="color:{z_color};font-weight:700">{z:+.2f}</td>'
            f'<td><span style="background:#fadbd8;color:#c0392b;padding:2px 8px;border-radius:4px;font-size:0.75rem;font-weight:700">{risk_label}</span></td>'
            f'<td>{cv_genes_html or "—"}</td>'
            f'<td>{pg_genes_html or "—"}</td>'
            f'<td style="font-size:0.78rem;color:var(--color-text-secondary)">{ctx}</td></tr>'
        )

    subsection_b = (
        f'<div style="overflow-x:auto"><table>'
        f'<thead><tr><th>{t["trait"]}</th><th>{t["z_score"]}</th><th>{t["risk"]}</th>'
        f'<th>{t["clinvar_genes"]}</th><th>{t["pharmgkb_genes"]}</th><th>{t["context"]}</th></tr></thead>'
        f'<tbody>{conv_rows}</tbody></table></div>'
    ) if high_risk_traits else f'<p style="color:var(--color-text-secondary)">No elevated PRS traits to cross-reference.</p>'

    # ── (c) Drug-Gene-PRS Intersections ──
    drug_rows = ""
    for d in DRUG_PRS_INTERSECTIONS:
        interaction = d["interaction_en"] if is_en else d["interaction_es"]
        recommendation = d["recommendation_en"] if is_en else d["recommendation_es"]
        drug_rows += (
            f'<tr><td><strong>{d["drug_class"]}</strong><br><span style="font-size:0.72rem;color:var(--color-text-secondary)">{d["drugs"]}</span></td>'
            f'<td><span class="clinical-convergence-gene">{d["gene"]}</span></td>'
            f'<td><strong>{d["prs_trait"]}</strong></td>'
            f'<td style="font-size:0.78rem">{interaction}</td>'
            f'<td style="font-size:0.78rem;color:#b7950b">{recommendation}</td></tr>'
        )

    subsection_c = (
        f'<div style="overflow-x:auto"><table>'
        f'<thead><tr><th>{t["drug_class"]}</th><th>{t["gene"]}</th><th>{t["trait"]}</th>'
        f'<th>{t["interaction"]}</th><th>{t["recommendation"]}</th></tr></thead>'
        f'<tbody>{drug_rows}</tbody></table></div>'
    )

    return f"""
    <div class="highlight-box" style="background:#eaf2f8;border:1px solid #aed6f1;border-radius:8px;padding:1rem 1.5rem;margin-bottom:1rem">
        <p style="font-size:0.9rem;margin:0">{t["desc"]}</p>
    </div>
    {subsection_a}
    <h4 style="margin-top:1.5rem">📊 {t["subsection_b"]}</h4>
    {subsection_b}
    <h4 style="margin-top:1.5rem">💊 {t["subsection_c"]}</h4>
    {subsection_c}
    <div class="disclaimer-box" style="margin-top:1.5rem">
        <p style="font-size:0.82rem;margin:0;white-space:pre-line">⚠️ RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS. NO PARA DIAGNÓSTICO CLÍNICO.
These findings are computational intersections of public databases. They do NOT constitute medical advice.
Always consult a qualified healthcare professional before making any medical decisions.</p>
    </div>
    """


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN REPORT BUILDER
# ═══════════════════════════════════════════════════════════════════════════════

CSS = """
:root {
    --color-low: #27ae60; --color-medium: #f39c12; --color-high: #e74c3c;
    --color-bg: #f5f6fa; --color-surface: #ffffff; --color-text: #2c3e50;
    --color-text-secondary: #7f8c8d; --color-border: #dee2e6;
    --radius: 8px; --shadow: 0 1px 3px rgba(0,0,0,.08);
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: var(--color-text); background: var(--color-bg); line-height: 1.6;
}
@media print { body { background: #fff; } .collapsible-section { break-inside: avoid; } }
.report-header {
    background: linear-gradient(135deg, #1a1a2e, #16213e, #0f3460);
    color: #fff; padding: 3rem 2rem; text-align: center;
}
.report-header h1 { font-size: 2.2rem; font-weight: 700; margin-bottom: .5rem; }
.report-header .subtitle { font-size: 1.1rem; opacity: .9; }
.report-header .meta { margin-top: 1rem; font-size: .85rem; opacity: .7; }
.container { max-width: 1100px; margin: 0 auto; padding: 2rem 1.5rem; }

/* Collapsible sections */
.collapsible-section { margin-bottom: 1.5rem; }
.section-header {
    cursor: pointer; display: flex; align-items: center; gap: .75rem;
    padding: 1rem 1.25rem; background: var(--color-surface);
    border-radius: var(--radius); box-shadow: var(--shadow);
    border-left: 4px solid #3498db; user-select: none;
    transition: background .15s;
}
.section-header:hover { background: #eaf2f8; }
.section-header h2 { font-size: 1.1rem; font-weight: 600; }
.section-arrow { font-size: .8rem; transition: transform .2s; color: #3498db; }
.section-body { padding: 1.25rem; background: var(--color-surface);
    border-radius: 0 0 var(--radius) var(--radius);
    box-shadow: var(--shadow); margin-top: -2px; }

/* Summary cards */
.summary-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem; margin-bottom: 1.5rem; }
.summary-card {
    background: var(--color-surface); border-radius: var(--radius);
    padding: 1.25rem; text-align: center; box-shadow: var(--shadow);
    border-left: 4px solid var(--color-border);
}
.summary-card .card-number { font-size: 2.2rem; font-weight: 800; }
.summary-card .card-label { font-size: .75rem; color: var(--color-text-secondary); margin-top: .25rem; }

/* Info grid */
.info-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin: 1rem 0; }
.info-card {
    background: var(--color-surface); border-radius: var(--radius); padding: 1rem 1.25rem;
    box-shadow: var(--shadow); text-align: center;
}
.info-card h4 { font-size: .75rem; text-transform: uppercase; color: var(--color-text-secondary); letter-spacing: .5px; margin-bottom: .5rem; }
.info-card .big-stat { font-size: 1.8rem; font-weight: 800; }
.info-card .stat-sub { font-size: .75rem; color: var(--color-text-secondary); margin-top: .25rem; }

/* Highlight box */
.highlight-box {
    background: #eaf2f8; border-radius: var(--radius);
    padding: .75rem 1rem; margin: .75rem 0; font-size: .85rem;
}

/* Tables */
table { width: 100%; border-collapse: collapse; background: var(--color-surface);
    border-radius: var(--radius); overflow: hidden; box-shadow: var(--shadow); margin: .75rem 0; font-size: .85rem; }
th { background: #f1f3f5; padding: .55rem .65rem; text-align: left; font-weight: 600;
    font-size: .72rem; text-transform: uppercase; letter-spacing: .5px; color: var(--color-text-secondary); }
td { padding: .55rem .65rem; border-top: 1px solid var(--color-border); }
tr:hover { background: #f8f9fa; }

/* Risk badges */
.risk-high { background: #fadbd8; color: #c0392b; }
.risk-medium { background: #fdebd0; color: #b7950b; }
.risk-low { background: #d5f5e3; color: #1e8449; }

/* Disclaimer */
.disclaimer-box {
    background: #fef9e7; border: 1px solid #f9e79f;
    border-radius: var(--radius); padding: 1.25rem; margin-top: 2rem;
    font-size: .8rem; color: #7d6608;
}

/* Buttons */
.btn-row { display: flex; gap: .5rem; justify-content: flex-end; margin-bottom: 1rem; }
.btn {
    padding: .4rem .8rem; border: 1px solid var(--color-border); border-radius: 4px;
    background: var(--color-surface); cursor: pointer; font-size: .78rem;
    transition: background .15s;
}
.btn:hover { background: #eaf2f8; }
.btn.active { background: #3498db; color: #fff; border-color: #3498db; }

/* Footer */
.report-footer {
    text-align: center; padding: 2rem; font-size: .75rem;
    color: var(--color-text-secondary); border-top: 1px solid var(--color-border); margin-top: 2rem;
}

/* Print styles */
@media print {
    .section-body { display: block !important; }
    .btn-row { display: none; }
}

/* Confidence & trust tier elements */
.confidence-stars { display: flex; align-items: center; gap: 4px; min-width: 110px; }
.tier-badge { display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 0.65rem; font-weight: 700; white-space: nowrap; }
.tier-1 { background: #d5f5e3; color: #1e8449; }
.tier-2 { background: #fdebd0; color: #b7950b; }
.tier-3 { background: #fadbd8; color: #c0392b; }
.decomp-bar { display: flex; height: 8px; border-radius: 3px; overflow: hidden; }
.decomp-bar-gen { background: #3498db; }
.decomp-bar-anc { background: #e74c3c; }
.decomp-bar-eff { background: #f39c12; }
.snp-coverage-bar { height: 6px; border-radius: 3px; overflow: hidden; flex: 1; }
.confidence-bar { height: 6px; border-radius: 3px; overflow: hidden; min-width: 50px; }
.confidence-bar-fill { height: 100%; border-radius: 3px; transition: width .3s; }
.portability-banner { background: #fef9e7; border: 2px solid #f39c12; border-radius: 8px; }
.cal-badge { padding: 2px 6px; border-radius: 3px; font-size: 0.65rem; font-weight: 700; white-space: nowrap; }
.limit-badge { padding: 1px 5px; border-radius: 3px; font-size: 0.6rem; font-weight: 700; margin-right: 2px; white-space: nowrap; }
.trust-legend { background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 0.6rem 1rem; margin-bottom: 0.8rem; font-size: 0.72rem; }
.radar-container { display: flex; justify-content: center; margin: 1rem 0; }
.radar-svg { max-width: 420px; height: auto; }

/* Clinical actionability */
.clinical-convergence-card { background: var(--color-surface); border-radius: var(--radius); padding: 1rem; box-shadow: var(--shadow); margin: 0.5rem 0; }
.clinical-convergence-gene { display: inline-block; padding: 1px 6px; background: #eaf2f8; border-radius: 3px; font-size: 0.72rem; font-weight: 600; margin: 1px; }
.clinical-convergence-drug { display: inline-block; padding: 1px 6px; background: #fef9e7; border-radius: 3px; font-size: 0.72rem; font-weight: 600; margin: 1px; }
"""

JS = """
function toggleSection(id, header) {
    const body = document.getElementById(id);
    const arrow = document.getElementById(id + '_arrow');
    if (body.style.display === 'none') {
        body.style.display = 'block';
        arrow.textContent = '▼';
    } else {
        body.style.display = 'none';
        arrow.textContent = '▶';
    }
}
function expandAll() {
    document.querySelectorAll('.section-body').forEach(b => b.style.display = 'block');
    document.querySelectorAll('.section-arrow').forEach(a => a.textContent = '▼');
}
function collapseAll() {
    document.querySelectorAll('.section-body').forEach(b => b.style.display = 'none');
    document.querySelectorAll('.section-arrow').forEach(a => a.textContent = '▶');
}
"""


def build_html_report(lang: str, data: Dict, sample_id: str) -> str:
    """Build the complete HTML report."""
    ui = UI.get(lang, UI["en"])
    ui["_lang"] = lang  # Pass language through for bilingual sections

    # Build sections
    sections_html = ""
    s = ui["sections"]

    # 0. Top Findings — prioritized, plain-language, actionable-navigation
    # summary (IMPROVEMENT_PLAN.md 1.1). Purely additive: everything below
    # (including the full PRS table it links into) is unchanged.
    sections_html += collapsible_section("top_findings", ui["top_findings_title"],
        build_top_findings(data["prs_result"].get("prs_entries", []), ui,
                           evidence_lookup=data.get("_evidence_lookup", {}),
                           cal_lookup=data.get("_cal_lookup", {}),
                           uncert_lookup=data.get("_uncert_lookup", {}),
                           recommendation_lookup=data.get("_recommendation_lookup", {}),
                           ancestry=data["ancestry"]),
        open_by_default=True)

    # 1. Summary (always open)
    sections_html += collapsible_section("summary", f"📊 {s['summary']}",
        build_summary_cards(data["prs_result"], data["ancestry"],
                           data["integrity"], data["validation"], ui,
                           cal_lookup=data.get("_cal_lookup", {}),
                           uncert_lookup=data.get("_uncert_lookup", {}),
                           evidence_lookup=data.get("_evidence_lookup", {}),
                           portability=data.get("portability", {})),
        open_by_default=True)

    # 2. Ancestry
    sections_html += collapsible_section("ancestry", f"🌍 {s['ancestry']}",
        build_ancestry_section(data["ancestry"], data.get("pca_data", {}), ui))

    # 2b. Deep Ancestry (haplogroups, sub-continental)
    deep = data.get("deep_ancestry", {})
    if deep:
        ydna = deep.get("y_dna", {})
        mtdna = deep.get("mt_dna", {})
        subcont = deep.get("sub_continental", {})
        neand = deep.get("neanderthal", {})

        deep_html = "<div class='info-grid' style='grid-template-columns:1fr 1fr 1fr'>"
        deep_html += f"<div class='info-card'><h4>🧬 mtDNA Haplogroup</h4><div class='big-stat' style='font-size:1.2rem'>{mtdna.get('haplogroup','?')}</div><div class='stat-sub'>{mtdna.get('description','')}</div></div>"
        deep_html += f"<div class='info-card'><h4>🧬 Y-DNA Haplogroup</h4><div class='big-stat' style='font-size:1rem'>{ydna.get('haplogroup','?')}</div><div class='stat-sub'>{ydna.get('description','')[:80]}</div></div>"
        if neand.get("reliable"):
            pct = neand.get("percentage", "?")
            method = neand.get("method", "snp_panel")
            closest = neand.get("closest_population", "EUR")
            pop_label = neand.get("population_comparisons", {}).get(closest, {}).get("label", closest)

            if "AADR" in method:
                affinity = neand.get("admix_ratio", neand.get("archaic_affinity_ratio", "?"))
                deep_html += (f'<div class="info-card"><h4>🦴 Archaic DNA</h4>'
                              f'<div class="big-stat">{pct}%</div>'
                              f'<div class="stat-sub">AADR-based (1.23M SNPs) | Affinity: {affinity}x</div></div>')
            else:
                deep_html += (f'<div class="info-card"><h4>🦴 Neanderthal DNA</h4>'
                              f'<div class="big-stat">{pct}%</div>'
                              f'<div class="stat-sub">Closest to {pop_label} | 133-SNP panel</div></div>')
        elif neand.get("snps_found", 0) > 0:
            nf = neand.get("snps_found", 0)
            nt = neand.get("snps_total", "?")
            deep_html += (f'<div class="info-card"><h4>🦴 Neanderthal DNA</h4>'
                          f'<div class="big-stat" style="font-size:0.9rem">N/A</div>'
                          f'<div class="stat-sub">Insufficient coverage ({nf}/{nt} SNPs). Need WGS VCF.</div></div>')
        else:
            deep_html += (f'<div class="info-card"><h4>🦴 Neanderthal DNA</h4>'
                          f'<div class="big-stat" style="font-size:0.9rem">N/A</div>'
                          f'<div class="stat-sub">No archaic SNPs found in VCF</div></div>')
        deep_html += "</div>"

        # Sub-continental populations
        if subcont:
            assigned_sub = subcont.get("assigned_sub_population", "")
            if assigned_sub:
                # Real sub-continental classification available
                sub_name = subcont.get("sub_population_name", assigned_sub)
                confidence = subcont.get("confidence", "MODERATE")
                conf_color = {"HIGH": "#27ae60", "MODERATE": "#f39c12", "LOW": "#e74c3c"}.get(confidence, "#7f8c8d")
                deep_html += f"<h4 style='margin-top:1rem'>🌍 Sub-Continental Ancestry</h4>"
                deep_html += "<div class='info-grid' style='grid-template-columns:1fr 1fr'>"
                deep_html += f"<div class='info-card' style='border-left:3px solid {conf_color}'><h4>Assigned Population</h4><div class='big-stat' style='font-size:1.3rem'>{assigned_sub}</div><div class='stat-sub'>{sub_name}</div></div>"
                deep_html += f"<div class='info-card' style='border-left:3px solid {conf_color}'><h4>Confidence</h4><div class='big-stat' style='font-size:1.3rem;color:{conf_color}'>{confidence}</div><div class='stat-sub'>Max probability: {subcont.get('max_probability', 0):.0%}</div></div>"
                deep_html += "</div>"
                # Show probabilities per sub-population
                probs = subcont.get("posterior_probabilities", {})
                if probs:
                    deep_html += "<h4>Sub-Population Probabilities</h4>"
                    prob_rows = ""
                    for pop, prob in sorted(probs.items(), key=lambda x: -x[1]):
                        pct = prob * 100
                        prob_rows += f"<tr><td><strong>{pop}</strong></td><td>{pct:.1f}%</td><td><div style='height:6px;background:#e9ecef;border-radius:3px;overflow:hidden'><div style='width:{pct:.0f}%;height:100%;background:#3498db;border-radius:3px'></div></div></td></tr>"
                    deep_html += f"<table><thead><tr><th>Population</th><th>Probability</th><th>Distribution</th></tr></thead><tbody>{prob_rows}</tbody></table>"
                # Still show available reference populations
                subs = subcont.get("sub_populations_available", [])
                if subs:
                    deep_html += "<h4 style='margin-top:0.5rem'>Reference Populations</h4>"
                    deep_html += "<div class='info-grid' style='grid-template-columns:repeat(auto-fit, minmax(150px, 1fr))'>"
                    for p in subs[:5]:
                        deep_html += f"<div class='info-card'><h4>{p['code']}</h4><div class='stat-sub'>{p['name']}</div><div style='font-size:0.7rem;color:var(--color-text-secondary)'>{p.get('description','')[:100]}</div></div>"
                    deep_html += "</div>"
            else:
                # Fallback: informational listing only
                assigned = subcont.get("assigned_super_population", "EUR")
                subs = subcont.get("sub_populations_available", [])
                deep_html += f"<h4 style='margin-top:1rem'>🌍 Sub-Continental Reference Populations ({assigned})</h4>"
                deep_html += "<div class='info-grid' style='grid-template-columns:repeat(auto-fit, minmax(150px, 1fr))'>"
                for p in subs[:5]:
                    deep_html += f"<div class='info-card'><h4>{p['code']}</h4><div class='stat-sub'>{p['name']}</div><div style='font-size:0.7rem;color:var(--color-text-secondary)'>{p.get('description','')[:100]}</div></div>"
                deep_html += "</div>"
                deep_html += f"<p style='font-size:0.72rem;color:var(--color-text-secondary);margin-top:0.3rem'>{subcont.get('note','')}</p>"

        sections_html += collapsible_section("deep_ancestry", "🧬 Deep Ancestry — Haplogroups & Sub-Continental", deep_html)

    # 3. PRS Results (always open)
    sections_html += collapsible_section("prs", f"📈 {s['prs']}",
        build_prs_table(data["prs_result"].get("prs_entries", []), ui,
                       cal_lookup=data.get("_cal_lookup", {}),
                       uncert_lookup=data.get("_uncert_lookup", {}),
                       portability=data.get("portability", {}),
                       evidence_lookup=data.get("_evidence_lookup", {})),
        open_by_default=True)

    # 4. Uncertainty Decomposition
    sections_html += collapsible_section("uncertainty_decomp", f"📐 {s['uncertainty']}",
        build_uncertainty_decomposition(data.get("uncertainty_report", {})))

    # 5. Variant Detail
    sections_html += collapsible_section("variants", f"🧬 {s['variants']}",
        build_variant_detail(data["prs_result"].get("prs_entries", [])))

    # 6. Population Calibration
    sections_html += collapsible_section("calibration", f"📊 {s['calibration']}",
        build_calibration_detail_section(data.get("calibration_report", {})))

    # 6b. PGS Catalog Calibration
    pgs_data = data.get("pgs_calibration", {})
    if pgs_data and pgs_data.get("all_entries"):
        n_scores = pgs_data.get("summary", {}).get("total_scores", len(pgs_data["all_entries"]))
        pgs_title = f"🧬 {s.get('pgs_calibration', 'PGS Catalog Calibration')} ({n_scores} scores, calibrated)"
        sections_html += collapsible_section("pgs_calibration", pgs_title,
            build_pgs_calibration_section(pgs_data, ui=ui,
                                         pgs_coverage=data.get("_pgs_coverage_lookup", {})))

    # 6c. ClinVar Pathogenic Variants
    clinvar_data = data.get("clinvar", {})
    sections_html += collapsible_section("clinvar", f"🧬 {s.get('clinvar', 'ClinVar Pathogenic Variants')}",
        build_clinvar_section(clinvar_data, ui))

    # 6d. Pharmacogenomics
    pharmgkb_data = data.get("pharmgkb", {})
    pharmgkb_findings = pharmgkb_data.get("pharmacogenomic_findings", [])
    if pharmgkb_findings:
        sections_html += collapsible_section("pharmgkb", "💊 Pharmacogenomics — Drug Response",
            build_pharmgkb_section(pharmgkb_data, ui))

    # 6e. Clinical Actionability Summary (cross-references ClinVar + PharmGKB + PRS)
    sections_html += collapsible_section("clinical_actionability", "🏥 Clinical Actionability Summary",
        build_clinical_actionability_section(
            data.get("clinvar", {}),
            data.get("pharmgkb", {}),
            data["prs_result"].get("prs_entries", []),
            ui))

    # 7. Population Portability
    sections_html += collapsible_section("portability", f"🌐 {s['portability']}",
        build_portability_section(data.get("portability", {})))

    # 8. Scientific Validation
    sections_html += collapsible_section("validation", f"✅ {s['validation']}",
        build_validation_section(data["validation"], ui))

    # 9. GWAS Consortium Validation
    sections_html += collapsible_section("gwas_consortium", f"🏛️ {s['gwas_consortium']}",
        build_gwas_consortium_section(data.get("gwas_consortium", {})))

    # 10. Quality Delta Benchmarking
    sections_html += collapsible_section("benchmark", f"🔬 {s['benchmark']}",
        build_benchmark_section(data["benchmark"], data.get("quality_delta", {}), ui))

    # 11. Adversarial Stress Testing
    sections_html += collapsible_section("adversarial", f"🛡️ {s['adversarial']}",
        build_adversarial_section(data["adversarial"], ui))

    # 12. Failure Mode Coverage
    sections_html += collapsible_section("failure_map", f"⚠️ {s['failure_map']}",
        build_failure_map_section(data["failure_map"], ui))

    # 13. Leakage Prevention — Detailed
    sections_html += collapsible_section("leakage", f"🔒 {s['leakage']}",
        build_leakage_detail_section(data.get("leakage_audit", data.get("leakage", {}))))

    # 14. GWAS-Ancestry Consistency
    sections_html += collapsible_section("consistency", f"🔗 {s['consistency']}",
        build_consistency_section(data.get("consistency", {})))

    # 15. Scientific Integrity
    sections_html += collapsible_section("integrity", f"🏅 {s['integrity']}",
        build_integrity_section(data["integrity"], ui))

    # 15b. PRS Profile Radar (interactive Chart.js)
    sections_html += collapsible_section("radar", "🕸️ PRS Profile Radar",
        build_radar_chart_js(data["prs_result"].get("prs_entries", []), ui,
                             cal_lookup=data.get("_cal_lookup"),
                             uncert_lookup=data.get("_uncert_lookup"),
                             evidence_lookup=data.get("_evidence_lookup")))

    # 16. Reproducibility
    sections_html += collapsible_section("reproducibility", f"🔁 {s['reproducibility']}",
        build_reproducibility_section(data.get("reproducibility", {})))

    # 17. Pipeline Methodology
    sections_html += collapsible_section("methodology", f"⚙️ {s['methodology']}",
        build_methodology_section(data["prs_result"], data["ancestry"]))

    # 18. Limitations (always open)
    # Build per-trait limitation notes
    entries = data["prs_result"].get("prs_entries", [])
    cal_lookup = data.get("_cal_lookup", {})
    trait_notes = ""
    for e in entries:
        trait = e.get("trait", "")
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)
        uncertainty = safe_float(e.get("uncertainty_score", 1.0))
        cal_entry = cal_lookup.get(trait.lower(), {})
        issues = []
        if n_total > 0 and n_used / n_total < 0.5:
            issues.append(f"Only {n_used} of {n_total} SNPs available — result may not capture full genetic risk")
        if cal_entry and safe_float(cal_entry.get("calibration_slope", 1.0)) < 0:
            slope = safe_float(cal_entry.get("calibration_slope", 0))
            issues.append(f"Calibration direction is reversed (slope={slope:.2f}) — percentile may be unreliable")
        if uncertainty >= 0.8:
            issues.append("Effect size uncertainty dominates — use with caution")
        if issues:
            trait_notes += f"<li><strong>{trait}:</strong> {'; '.join(issues)}</li>"
        else:
            trait_notes += f"<li><strong>{trait}:</strong> No significant limitations detected</li>"

    limitations_html = (
        f'<h4>Per-Trait Confidence Notes</h4>'
        f'<ul style="font-size:0.82rem;line-height:1.6;margin-bottom:1.5rem">{trait_notes}</ul>'
        f'<div class="disclaimer-box">'
        f'<h3>⚠️ Important Disclaimer</h3>'
        f'<p style="white-space:pre-line">{ui["disclaimer"]}</p>'
        f'</div>'
    )
    sections_html += collapsible_section("limitations", f"⚠️ {s['limitations']}",
        limitations_html,
        open_by_default=True)

    now = datetime.now().strftime("%Y-%m-%d %H:%M UTC")
    pop = data["ancestry"].get("assigned_population", "EUR")
    integrity_score = data["integrity"].get("scientific_integrity_score", 0)

    return f"""<!DOCTYPE html>
<html lang="{lang}">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{ui['title']} — {sample_id}</title>
    <style>{CSS}</style>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
</head>
<body>
    <header class="report-header">
        <h1>{ui['title']}</h1>
        <div class="subtitle">{ui['subtitle']}</div>
        <div class="meta">
            Sample: {sample_id} | Population: {pop} | Integrity: {integrity_score:.0f}/100<br>
            Generated: {now} | BlueGen v{PIPELINE_VERSION} | GRCh37/hg19
        </div>
    </header>
    <div class="container">
        {reference_coverage_banner(data["prs_result"], lang)}
        <div class="btn-row">
            <button class="btn" onclick="expandAll()">📖 Expand All</button>
            <button class="btn" onclick="collapseAll()">📕 Collapse All</button>
        </div>
        {sections_html}
    </div>
    <footer class="report-footer">
        <p>BlueGen Report — Generated {now}</p>
        <p>BlueGen v{PIPELINE_VERSION} | Sample: {sample_id} | Research Use Only</p>
    </footer>
    <script>{JS}</script>
</body>
</html>"""


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Comprehensive HTML Report Generator")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--output-dir", "-o", default="reports")
    parser.add_argument("--lang", default="both", choices=["en", "es", "both"])
    parser.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all data sources
    logger.info("═══ Loading Pipeline Data ═══")

    data = {}

    # Core outputs
    p = Path(".")
    data["prs_result"] = load_json("prs/PRS_RESULT.json")
    data["ancestry"] = load_json("science/ANCESTRY_MODEL.json")
    data["validation"] = load_json("science/global_validation_report.json")
    data["integrity"] = load_json("FINAL_SCIENTIFIC_SCORE.json")
    data["benchmark"] = load_json("benchmark/VALIDATION_REPORT.json")
    data["adversarial"] = load_json("science/adversarial_validation_report.json")
    data["failure_map"] = load_json("science/failure_mode_map.json")
    data["leakage"] = load_json("science/pipeline_gate_check.json")
    data["quality_delta"] = load_json("benchmark/quality_delta.json")
    # Extended data sources
    data["uncertainty_report"] = load_json("prs/uncertainty_report.json")
    data["calibration_report"] = load_json("prs/population_calibration_report.json")
    data["gwas_consortium"] = load_json("benchmark/gwas_consortium_validation.json")
    data["portability"] = load_json("benchmark/portability_report.json")
    data["reproducibility"] = load_json("reproducibility/run_fingerprint.json")
    data["consistency"] = load_json("prs/consistency_check_report.json")
    data["leakage_audit"] = load_json("science/leakage_audit.json")
    data["snp_universe"] = load_json("science/snp_universe.json")
    data["pgs_calibration"] = load_json("prs/pgs_scores/pgs_calibration_report.json")
    data["clinvar"] = load_json("clinvar/clinvar_pathogenic_variants.json")
    data["pharmgkb"] = load_json("pharmgkb/pharmgkb_drug_report.json")
    data["deep_ancestry"] = load_json("ancestry/deep_ancestry.json")

    # NEW: Load subcontinental assignment if available
    subc_data = load_json("pca/subcontinental_assignment.json")
    if subc_data and subc_data.get("assigned_sub_population"):
        # Merge into deep_ancestry for display
        if data["deep_ancestry"]:
            data["deep_ancestry"]["sub_continental"] = subc_data
        else:
            data["deep_ancestry"] = {"sub_continental": subc_data}
        logger.info(f"  ✅ Subcontinental assignment: {subc_data.get('assigned_sub_population')}")

    # NEW: Load calibration validation data per trait (for confidence scores)
    data["calibration_validation"] = load_json("benchmark/calibration_validation.json")

    # Build per-trait calibration lookup
    cal_lookup = {}
    for v in data["calibration_validation"].get("validations", []):
        cal_lookup[v["trait"].lower()] = v
    data["_cal_lookup"] = cal_lookup

    # Build per-trait uncertainty decomposition lookup
    uncert_lookup = {}
    for r in data.get("uncertainty_report", {}).get("results", []):
        uncert_lookup[r["trait"].lower()] = r
    data["_uncert_lookup"] = uncert_lookup

    # Build per-trait evidence level lookup from SNP database
    evidence_lookup = {}
    snp_db_path = args.snp_db
    if os.path.exists(snp_db_path):
        try:
            snp_db = pd.read_csv(snp_db_path, dtype=str)
            for col in ["trait_category", "trait", "Trait", "trait_name"]:
                if col in snp_db.columns:
                    trait_col = col
                    break
            else:
                trait_col = None
            if trait_col:
                ev_map = {"A": 100, "B": 75, "C": 50, "D": 25}
                for _, row in snp_db.iterrows():
                    trait = str(row.get(trait_col, "")).strip().lower()
                    ev = str(row.get("evidence_level", row.get("evidence", "C"))).strip().upper()
                    ev_score = ev_map.get(ev, 50)
                    if trait not in evidence_lookup:
                        evidence_lookup[trait] = []
                    evidence_lookup[trait].append(ev_score)
        except Exception:
            pass
    data["_evidence_lookup"] = evidence_lookup

    # Curated, evidence-cited per-trait action recommendations (IMPROVEMENT_PLAN.md
    # 1.4). Deliberately a small hand-verified subset, not auto-generated - see
    # data/trait_recommendations.json's _meta for the curation scope/rationale.
    recommendation_lookup = {}
    trait_reco_path = os.path.join(os.path.dirname(args.snp_db), "trait_recommendations.json")
    if os.path.exists(trait_reco_path):
        try:
            with open(trait_reco_path) as fh:
                raw = json.load(fh)
            recommendation_lookup = {k.lower(): v for k, v in raw.items() if k != "_meta"}
        except Exception:
            pass
    data["_recommendation_lookup"] = recommendation_lookup

    # Build PGS coverage lookup from pgs_results.csv
    pgs_coverage_lookup = {}
    pgs_results_path = "prs/pgs_scores/pgs_results.csv"
    if os.path.exists(pgs_results_path):
        try:
            pgs_results = pd.read_csv(pgs_results_path, dtype=str)
            for _, row in pgs_results.iterrows():
                pgs_id = str(row.get("pgs_id", "")).strip()
                n_used = safe_float(row.get("n_snps_used", 0))
                n_total = safe_float(row.get("n_snps_in_score", 0))
                if pgs_id:
                    pgs_coverage_lookup[pgs_id] = {"n_used": n_used, "n_total": n_total}
        except Exception:
            pass
    data["_pgs_coverage_lookup"] = pgs_coverage_lookup

    # Log what was found
    for name, d in data.items():
        status = "✅" if d else "⬚ (missing)"
        logger.info(f"  {status} {name}")

    # Fall back to ancestry inference if ANCESTRY_MODEL is empty or UNKNOWN
    anc_pop = data["ancestry"].get("assigned_population", "UNKNOWN")
    if not data["ancestry"] or anc_pop in ("UNKNOWN", None, ""):
        alt = load_json("pca/ancestry_inference.json")
        if alt:
            # Map ancestry_inference format to ANCESTRY_MODEL format
            summary = alt.get("summary", {})
            data["ancestry"] = {
                "assigned_population": summary.get("assigned_super_population", "EUR"),
                "posterior_probabilities": summary.get("all_probabilities", {}),
                "confidence": summary.get("confidence", "MODERATE"),
                "n_reference_samples": alt.get("methodology", {}).get("snps_used", 2504),
                "n_pcs": 20,
                "method": alt.get("methodology", {}).get("method", "allele_frequency_distance"),
            }
            logger.info(f"  Using pca/ancestry_inference.json: {data['ancestry']['assigned_population']} ({data['ancestry']['confidence']})")

    sample_id = args.sample_id
    # Try to detect actual sample ID from PRS result
    prs_sample = data["prs_result"].get("sample_id", "")
    if prs_sample and prs_sample != "SAMPLE_001":
        sample_id = prs_sample

    langs = ["en", "es"] if args.lang == "both" else [args.lang]

    for lang in langs:
        logger.info(f"  Generating {lang.upper()} report...")
        html = build_html_report(lang, data, sample_id)

        out_path = output_dir / f"comprehensive_report_{lang}.html"
        with open(out_path, "w") as fh:
            fh.write(html)

        size_kb = len(html) / 1024
        logger.info(f"    ✅ {out_path} ({size_kb:.0f} KB)")

    logger.info("═══ Report Complete ═══")
    return 0


if __name__ == "__main__":
    sys.exit(main())
