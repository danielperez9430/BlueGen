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
from clinical.disease_taxonomy import BODY_SYSTEM_ORDER, system_label

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
        # For traits where "high" genuinely means a favorable direction (e.g. Morning
        # chronotype, Cognitive function - see _polarity_inverted in trait_recommendations.json)
        # rather than elevated health risk. Never uses red/"RISK" wording, since these
        # traits don't represent danger even at their "low" end.
        "favorable_high": "FAVORABLE", "favorable_medium": "TYPICAL", "favorable_low": "LESS FAVORABLE",
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
        "favorable_high": "FAVORABLE", "favorable_medium": "TÍPICO", "favorable_low": "MENOS FAVORABLE",
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
# DATA LOADING + INTERPRETATION HELPERS
# (extracted to report/interpretations.py and report/data_loader.py -
#  IMPROVEMENT_PLAN.md 1.6 Phase 1 - imported here verbatim, no logic changes)
# ═══════════════════════════════════════════════════════════════════════════════

from report.interpretations import (
    load_json, safe_float, trait_anchor_id, risk_color, risk_badge, risk_bar,
    compute_per_trait_confidence, confidence_stars, calibration_flag, trust_tier,
    trust_badge, mini_decomp_bar, snp_coverage_bar, portability_banner,
    reference_coverage_banner, trait_limitations_badges, trust_tier_legend,
    evidence_letter, evidence_badge,
)
from report.data_loader import load_report_data
from report.render import build_html_report as render_document, render_partial

# ═══════════════════════════════════════════════════════════════════════════════
# RADAR CHART (Chart.js — interactive)
# ═══════════════════════════════════════════════════════════════════════════════

def build_radar_chart_js(entries, ui, cal_lookup=None, uncert_lookup=None, evidence_lookup=None, polarity_inverted=None):
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
    if polarity_inverted is None:
        polarity_inverted = set()

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
        color = risk_color(z, inverted=trait_name.lower() in polarity_inverted)
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

def build_top_findings(entries, ui, evidence_lookup=None, cal_lookup=None, uncert_lookup=None,
                        recommendation_lookup=None, ancestry=None, max_findings=8,
                        polarity_inverted=None):
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
    if polarity_inverted is None: polarity_inverted = set()
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
        return render_partial("top_findings.html.j2", empty_message=ui["top_findings_empty"])

    risk_words = {"high": ui["risk_word_high"], "medium": ui["risk_word_medium"], "low": ui["risk_word_low"]}
    cards = []
    for priority, e, z, pctl, risk, ev_avg, conf, n_used, n_total in top:
        trait = e.get("trait", "")
        meaning = ui["top_findings_meaning"].format(
            trait=trait, risk_word=risk_words.get(risk, risk_words["medium"]),
            pctl=pctl, pop=pop, z=z, n_used=n_used, n_total=n_total,
            evidence=evidence_letter(ev_avg))
        lang = ui.get("_lang", "en")
        curated = recommendation_lookup.get(trait.lower()) if risk == "high" else None
        recommendation = (curated or {}).get("recommendation_" + lang) or ui["top_findings_action_fallback"]
        inverted = trait.lower() in polarity_inverted
        cards.append({
            "color": risk_color(z, inverted=inverted), "trait": trait,
            "risk_badge": risk_badge(risk, ui, inverted=inverted), "evidence_badge": evidence_badge(ev_avg),
            "meaning": meaning, "recommendation": recommendation, "anchor": trait_anchor_id(trait),
        })

    # Compact, collapsible confidence-context note (IMPROVEMENT_PLAN.md 1.5)
    # — the full disclaimer + per-trait confidence notes already exist at
    # the bottom of the report (Limitations & Disclaimers section); this is
    # the same text, just also reachable without scrolling past everything.
    return render_partial("top_findings.html.j2",
        empty_message=None, intro=ui["top_findings_intro"],
        disclaimer_summary=ui["top_findings_disclaimer_summary"], disclaimer_text=ui["disclaimer"],
        cards=cards, jump_label=ui["top_findings_jump"])

def build_summary_cards(prs_result, ancestry, integrity, validation, ui, cal_lookup=None, uncert_lookup=None, evidence_lookup=None, portability=None, polarity_inverted=None):
    """Executive summary cards with confidence overview."""
    if cal_lookup is None:
        cal_lookup = {}
    if uncert_lookup is None:
        uncert_lookup = {}
    if evidence_lookup is None:
        evidence_lookup = {}
    if polarity_inverted is None:
        polarity_inverted = set()

    entries = prs_result.get("prs_entries", [])
    # For polarity-inverted traits (risk_category=="high" means a favorable
    # result, not elevated risk - see risk_badge docstring), count them under
    # the opposite bucket so this "n traits at higher/lower risk" summary
    # isn't itself misleading. "medium" is unaffected - polarity only flips
    # which end is favorable, not the middle.
    def _bucket(e):
        risk = e.get("risk_category")
        if risk in ("high", "low") and e.get("trait", "").lower() in polarity_inverted:
            return "low" if risk == "high" else "high"
        return risk
    n_high = sum(1 for e in entries if _bucket(e) == "high")
    n_medium = sum(1 for e in entries if _bucket(e) == "medium")
    n_low = sum(1 for e in entries if _bucket(e) == "low")

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

    top_trait_ctx = None
    if top_trait:
        top_trait_ctx = {
            "trait": top_trait["trait"], "raw_score": f"{safe_float(top_trait.get('raw_score', 0)):.2f}",
            "n_used": top_trait.get("n_snps_used", 0), "n_total": top_trait.get("n_snps_total", 0),
        }

    return render_partial("summary_cards.html.j2",
        n_high=n_high, risk_high_label=ui['risk_high'],
        n_medium=n_medium, risk_medium_label=ui['risk_medium'],
        n_low=n_low, risk_low_label=ui['risk_low'],
        conf_color=conf_color, avg_conf=f"{avg_conf:.0f}",
        tier1_count=tier_counts.get('TIER 1', 0), tier3_count=tier_counts.get('TIER 3', 0),
        pop_name=pop_name, confidence=confidence, integrity_score=f"{integrity_score:.0f}",
        top_trait=top_trait_ctx, integrity_cat=integrity_cat,
        integrity_cat_desc=integrity.get("category_description", ""),
        best_trait=best_trait, best_conf=f"{best_conf:.0f}",
        worst_trait=worst_trait, worst_conf=f"{worst_conf:.0f}",
        port_note=port_note)


def build_ancestry_section(ancestry, pca_data, ui):
    """Ancestry deep-dive section."""
    pop = ancestry.get("assigned_population", "UNKNOWN")
    confidence = ancestry.get("confidence", "UNKNOWN")
    probs = ancestry.get("posterior_probabilities", {})
    n_ref = ancestry.get("n_reference_samples", 2504)
    n_pcs = ancestry.get("n_pcs", 20)

    pop_names = POP_NAMES["en"]
    prob_rows = []
    for p in ["EUR", "AFR", "EAS", "SAS", "AMR"]:
        prob = probs.get(p, 0) * 100 if isinstance(probs.get(p, 0), (int, float)) else 0
        prob_rows.append({
            "label": f"{pop_names.get(p, p)} ({p})",
            "pct": f"{prob:.1f}",
            "bar_w": max(prob, 1),
            "bar_color": "#3498db" if p == pop else "#bdc3c7",
        })

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

    return render_partial("ancestry.html.j2",
        pop_display=pop_names.get(pop, pop), confidence=confidence,
        n_ref=n_ref, n_pcs=n_pcs, prob_rows=prob_rows, pca_table_html=pca_table)


def build_prs_table(entries, ui, cal_lookup=None, uncert_lookup=None, portability=None, evidence_lookup=None,
                     polarity_inverted=None):
    """Full PRS results table with risk bars and confidence metrics."""
    if cal_lookup is None:
        cal_lookup = {}
    if uncert_lookup is None:
        uncert_lookup = {}
    if evidence_lookup is None:
        evidence_lookup = {}
    if polarity_inverted is None:
        polarity_inverted = set()

    rows = []
    confidences = []
    tier_counts = {"TIER 1": 0, "TIER 2": 0, "TIER 3": 0}

    for e in entries:
        trait = e.get("trait", "")
        z = safe_float(e.get("population_zscore", e.get("raw_score", 0)))
        pctl = safe_float(e.get("population_percentile", 50))
        risk = e.get("risk_category", "medium")
        n_used = e.get("n_snps_used", 0)
        n_total = e.get("n_snps_total", 0)
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
        inverted = trait.lower() in polarity_inverted
        color = risk_color(z, inverted=inverted)

        rows.append({
            "anchor_id": trait_anchor_id(trait), "trait": trait, "color": color, "z": f"{z:+.2f}",
            "pctl": f"{pctl:.1f}", "risk_badge": risk_badge(risk, ui, inverted=inverted),
            "confidence_stars": confidence_stars(conf_score), "calibration_flag": calibration_flag(cal_entry),
            "trust_badge": trust_badge(tier), "snp_coverage_bar": snp_coverage_bar(n_used, n_total),
            "mini_decomp_bar": mini_decomp_bar(decomp),
            "limitations_badges": trait_limitations_badges(e, cal_entry),
            "risk_bar": risk_bar(bar_pct, z, inverted=inverted),
        })

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

    return render_partial("prs_table.html.j2",
        port_banner=port_banner, trust_tier_legend_html=trust_tier_legend(),
        summary_bar=summary_bar, rows=rows)


def build_variant_detail(entries, snp_db_path="data/snp_database_annotated.csv"):
    """Per-trait variant-level detail tables."""
    # Try to load the SNP database
    snp_db = None
    if os.path.exists(snp_db_path):
        try: snp_db = pd.read_csv(snp_db_path, dtype=str)
        except Exception: pass

    ev_colors = {"A": ("#27ae60", "#d5f5e3"), "B": ("#2e86c1", "#d6eaf8"), "C": ("#f39c12", "#fdebd0"), "D": ("#95a5a6", "#eaecee")}
    trait_sections = []
    for e in entries:
        trait = e.get("trait", "")

        # Detect trait column name (snake_case variants)
        trait_col = None
        if snp_db is not None:
            for col in ["trait_category", "trait", "Trait", "trait_name"]:
                if col in snp_db.columns:
                    trait_col = col
                    break

        if trait_col:
            trait_snps = snp_db[snp_db[trait_col].str.lower() == trait.lower()]
        else:
            trait_snps = None

        variants = []
        if trait_snps is not None and len(trait_snps) > 0:
            for _, row in trait_snps.iterrows():
                evidence = row.get("evidence_level", row.get("evidence", "—")).strip().upper()
                fg, bg = ev_colors.get(evidence, ("#7f8c8d", "#e9ecef"))
                variants.append({
                    "rsid": row.get("rsid", "—"), "gene": row.get("gene", "—"),
                    "effect_allele": row.get("effect_allele", "—"),
                    "weight": row.get("weight", row.get("beta", "—")),
                    "ev_badge": f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.7rem;font-weight:700">{evidence}</span>',
                })

        trait_sections.append({
            "trait": trait, "n_used": e.get("n_snps_used", 0), "n_total": e.get("n_snps_total", 0),
            "variants": variants,
        })

    return render_partial("variant_detail.html.j2", trait_sections=trait_sections)


def build_validation_section(validation, ui):
    """Validation checks table."""
    checks = validation.get("checks", [])
    sev_colors = {"ERROR": "#e74c3c", "WARNING": "#f39c12", "INFO": "#3498db"}
    rows = [{
        "check_id": c.get("check_id", ""), "category": c.get("category", ""),
        "description": c.get("description", ""),
        "sev": c.get("severity", "INFO"),
        "sev_color": sev_colors.get(c.get("severity", "INFO"), "#7f8c8d"),
        "icon": "✅" if c.get("passed", False) else "❌",
        "detail": c.get("detail", ""),
    } for c in checks]

    return render_partial("validation.html.j2",
        overall_score=f"{validation.get('overall_score', 0):.0f}",
        overall_status=validation.get('overall_status', '').replace('_', ' ').title(),
        passed=validation.get('passed', 0), total_checks=validation.get('total_checks', 0),
        warnings=validation.get('warnings', 0), errors=validation.get('errors', 0),
        rows=rows)


def build_benchmark_section(benchmark, quality_delta, ui):
    """GWAS & external benchmarking."""
    entries = benchmark.get("entries", [])
    bench_rows = []
    for e in entries:
        is_circ = e.get("is_circular", False)
        circ_badge = '<span style="background:#fdebd0;color:#b7950b;padding:1px 5px;border-radius:3px;font-size:0.65rem">CIRCULAR</span>' if is_circ else '<span style="background:#d5f5e3;color:#1e8449;padding:1px 5px;border-radius:3px;font-size:0.65rem">INDEPENDENT</span>'
        bench_rows.append({
            "validation_id": e.get('validation_id', ''), "description": e.get('description', ''),
            "vtype": e.get("validation_type", "unknown"), "circ_badge": circ_badge,
            "status": e.get("status", "VALID"),
        })

    # Quality delta
    qd = quality_delta
    components = qd.get("components", [])
    delta_rows = []
    for c in components:
        delta = c.get("delta", 0)
        direction = c.get("direction", "at_par")
        d_color = "#27ae60" if direction == "overperform" else ("#e74c3c" if direction == "underperform" else "#f39c12")
        delta_rows.append({
            "dimension": c.get('dimension', ''),
            "internal_score": f"{safe_float(c.get('internal_score', 0)):.0f}",
            "external_benchmark": f"{safe_float(c.get('external_benchmark', 0)):.0f}",
            "d_color": d_color, "delta": f"{delta:+.0f}",
            "direction": direction.replace('_', ' ').title(),
            "explanation": c.get('explanation', '')[:150],
        })

    return render_partial("benchmark.html.j2",
        vs=benchmark.get('validation_summary', {}), bench_rows=bench_rows,
        mean_delta_color="#27ae60" if qd.get('mean_delta', 0) >= 0 else "#e74c3c",
        mean_delta=f"{qd.get('mean_delta', 0):+.1f}", qd=qd, delta_rows=delta_rows)


def build_adversarial_section(adversarial, ui):
    """Adversarial stress testing results."""
    results = adversarial.get("results", [])
    sev_color_map = {"CRITICAL": "#e74c3c", "HIGH": "#e67e22", "MODERATE": "#f39c12"}
    rows = []
    for r in results:
        robust = r.get("is_robust", False)
        if isinstance(robust, str):
            robust = robust.lower() == "true"
        sev = r.get("severity", "MODERATE")
        rows.append({
            "test_id": r.get('test_id', ''), "description": r.get('description', ''),
            "sev": sev, "sev_color": sev_color_map.get(sev, "#7f8c8d"),
            "icon": "✅" if robust else "❌",
            "robust_label": "Robust" if robust else "Vulnerable",
            "change": f"{safe_float(r.get('relative_change', 0)):+.2f}",
            "detail": r.get('detail', ''),
        })

    critical_findings = adversarial.get("critical_findings", [])

    return render_partial("adversarial.html.j2",
        score=f"{adversarial.get('overall_robustness_score', 0):.0f}",
        n_tests=adversarial.get('n_tests', 0), n_robust=adversarial.get('n_robust', 0),
        n_vulnerable=adversarial.get('n_vulnerable', 0),
        critical_findings=", ".join(critical_findings) if critical_findings else "",
        rows=rows)


def build_failure_map_section(failure_map, ui):
    """Failure mode coverage."""
    failures = failure_map.get("failures", [])
    sev_color_map = {"CRITICAL": "#e74c3c", "HIGH": "#e67e22", "MODERATE": "#f39c12"}
    rows = []
    for f in failures:
        sev = f.get("severity", "MODERATE")
        validated = f.get("adversarial_validated", False)
        rows.append({
            "id": f.get('id', ''), "component": f.get('component', ''),
            "failure": f.get('failure', ''), "sev": sev,
            "sev_color": sev_color_map.get(sev, "#7f8c8d"),
            "v_icon": "✅" if validated else "⬚",
            "effect": f.get('effect', '')[:120],
        })

    return render_partial("failure_map.html.j2",
        n_failures=failure_map.get('n_failures', 0), n_critical=failure_map.get('n_critical', 0),
        n_high=failure_map.get('n_high', 0),
        most_vulnerable=failure_map.get('most_vulnerable_component', 'N/A'), rows=rows)


def build_integrity_section(integrity, ui):
    """Scientific integrity score breakdown."""
    components = integrity.get("components", [])
    rows = []
    for c in components:
        score = safe_float(c.get("score", 0))
        rows.append({
            "name": c.get('name', ''),
            "score": f"{score:.1f}",
            "score_color": "#27ae60" if score >= 80 else ("#f39c12" if score >= 60 else "#e74c3c"),
            "weight": f"{safe_float(c.get('weight', 0))*100:.0f}",
            "contribution": f"{safe_float(c.get('contribution', 0)):.1f}",
            "source": c.get('source', ''),
        })

    total = integrity.get("scientific_integrity_score", 0)
    cat = integrity.get("category", "Unknown")
    cat_color = {"PUBLICATION_READY": "#27ae60", "RESEARCH_GRADE": "#3498db",
                 "NEEDS_REVISION": "#f39c12", "SIGNIFICANT_ISSUES": "#e67e22",
                 "NOT_PUBLISHABLE": "#e74c3c"}.get(cat, "#7f8c8d")

    return render_partial("integrity.html.j2",
        total=f"{total:.1f}", cat_color=cat_color, cat_display=cat.replace('_', ' ').title(),
        formula=integrity.get('formula', ''),
        weights_locked_icon="✅" if integrity.get('weights_locked', False) else "❌", rows=rows)


def build_uncertainty_decomposition(uncertainty_report):
    """Per-trait variance decomposition: genotype vs ancestry vs effect."""
    results = uncertainty_report.get("results", [])
    if not results:
        return "<p style='color:#7f8c8d'>Uncertainty report not available.</p>"

    rows = []
    for r in results:
        decomp = r.get("decomposition", {})
        gen_frac = safe_float(decomp.get("genotype_fraction", 0)) * 100
        anc_frac = safe_float(decomp.get("ancestry_fraction", 0)) * 100
        eff_frac = safe_float(decomp.get("effect_fraction", 0)) * 100
        rows.append({
            "trait": r.get("trait", ""),
            "prs": f"{safe_float(r.get('prs_point_estimate', 0)):.3f}",
            "se": f"{safe_float(r.get('prs_std_error', 0)):.3f}",
            "total_var": f"{safe_float(decomp.get('total_variance', 0)):.4f}",
            "gen_frac": f"{gen_frac:.1f}", "anc_frac": f"{anc_frac:.1f}", "eff_frac": f"{eff_frac:.1f}",
            "gen_frac_int": f"{gen_frac:.0f}", "anc_frac_int": f"{anc_frac:.0f}", "eff_frac_int": f"{eff_frac:.0f}",
            "n_genotype": r.get('n_snps_with_genotype', 0), "n_effect_se": r.get('n_snps_with_effect_se', 0),
        })

    return render_partial("uncertainty_decomposition.html.j2", rows=rows)


def build_gwas_consortium_section(gwas_consortium):
    """GWAS consortium validation — consortia info + trait-level checks."""
    consortia = gwas_consortium.get("consortia", {})
    validations = gwas_consortium.get("validations", [])

    cons_cards = [{
        "name": name, "primary_ancestry": c.get('primary_ancestry', 'EUR'),
        "n_discovery": f"{c.get('n_discovery', 0):,}", "pmid": c.get('pmid', ''),
        "traits": ', '.join(c.get('traits', [])[:3]),
    } for name, c in consortia.items()]

    val_rows = []
    for v in validations:
        passed = v.get("overall_status") == "PASS"
        val_rows.append({
            "consortium": v.get('consortium', ''), "trait": v.get('trait', ''),
            "match": f"{safe_float(v.get('effect_direction_match', 0)) * 100:.0f}",
            "overlap_count": v.get('snp_overlap_count', 0),
            "overlap_pct": f"{safe_float(v.get('snp_overlap_pct', 0)) * 100:.1f}",
            "icon": "✅" if passed else "❌", "status": v.get('overall_status', ''),
        })

    return render_partial("gwas_consortium.html.j2",
        n_consortia=len(consortia), total_count=gwas_consortium.get("total_checks", 0),
        passed_count=gwas_consortium.get("passed", 0), failed_count=gwas_consortium.get('failed', 0),
        cons_cards=cons_cards, val_rows=val_rows)


def build_portability_section(portability):
    """Population portability — PRS shift across populations."""
    pops = portability.get("populations", [])
    status_color_map = {"GOOD_PORTABILITY": "#27ae60", "MODERATE_PORTABILITY": "#f39c12", "LIMITED_PORTABILITY": "#e74c3c"}
    rows = []
    for p in pops:
        status = p.get("status", "")
        rows.append({
            "population": p.get("population", ""), "n_ref": p.get('n_reference_samples', 0),
            "prs_shift": f"{safe_float(p.get('mean_prs_shift', 0)):.2f}",
            "calib_drift": f"{safe_float(p.get('calibration_drift', 0)):.2f}",
            "rank_instability": f"{safe_float(p.get('rank_instability', 0)):.2f}",
            "bias_index": f"{safe_float(p.get('ancestry_bias_index', 0)):.3f}",
            "color": status_color_map.get(status, "#7f8c8d"),
            "status": status.replace('_', ' ').title(),
        })

    return render_partial("portability.html.j2",
        global_bias=f"{safe_float(portability.get('global_bias_index', 0)):.3f}",
        most_biased=portability.get('most_biased', 'N/A'),
        least_biased=portability.get('least_biased', 'N/A'), rows=rows)


def build_reproducibility_section(repro):
    """Reproducibility: environment, seeds, versions."""
    env = repro.get("environment", {})
    seeds = repro.get("seeds", {})

    tools = [{"tool": t, "version": v} for t, v in env.get("system_tools", {}).items()]
    packages = [{"pkg": pkg, "ver": ver} for pkg, ver in sorted(env.get("pip_packages", {}).items())[:20]]

    return render_partial("reproducibility.html.j2",
        run_id=repro.get('run_id', 'N/A')[:16],
        repro_score=f"{repro.get('reproducibility_score', 0):.0f}",
        pipeline_version=repro.get('pipeline_version', 'N/A'),
        global_seed=seeds.get('global_seed', 'N/A'),
        global_seed_table=seeds.get('global_seed', ''),
        os_name=env.get('os_name', ''), os_version=env.get('os_version', ''),
        architecture=env.get('architecture', ''),
        python_version=env.get('python_version', ''),
        python_implementation=env.get('python_implementation', ''),
        kernel=env.get('kernel', '')[:120],
        tools=tools, packages=packages,
        numpy_seed=seeds.get('numpy_seed', ''), python_hash_seed=seeds.get('python_hash_seed', ''),
        sklearn_seed=seeds.get('sklearn_seed', ''), plink_seed=seeds.get('plink_seed', ''),
        bootstrap_seed_0=seeds.get('bootstrap_seeds', [0])[0])


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
    rows = [{
        "trait": t.get('trait', ''), "gwas_pop": t.get("gwas_population", ""),
        "gwas_type": t.get("gwas_type", "").replace('_', ' ').title(),
        "target_pop": t.get('target_population', ''),
        "icon": "✅" if t.get("is_match", False) else "❌", "note": t.get('note', ''),
    } for t in trait_checks]

    return render_partial("consistency.html.j2",
        gwas_match_color="#27ae60" if consistency.get('gwas_ancestry_match') else "#e74c3c",
        gwas_match_label="✅ PASS" if consistency.get('gwas_ancestry_match') else "❌ FAIL",
        ld_match_color="#27ae60" if consistency.get('ld_ancestry_match') else "#f39c12",
        ld_match_label="✅ PASS" if consistency.get('ld_ancestry_match') else "⚠️ WARN",
        confidence_downgrade=f"{safe_float(consistency.get('confidence_downgrade', 0)):.2f}",
        recommended_gwas=consistency.get('recommended_gwas_source', 'N/A'),
        gwas_summary_html=_build_gwas_summary(trait_checks), rows=rows)


def build_leakage_detail_section(leakage_audit):
    """Detailed leakage audit — 7 checks."""
    checks = leakage_audit.get("checks", [])
    sev_color_map = {"ERROR": "#e74c3c", "WARNING": "#f39c12", "INFO": "#3498db"}
    rows = [{
        "check_id": c.get('check_id', ''), "description": c.get('description', ''),
        "sev": c.get("severity", "INFO"),
        "sev_color": sev_color_map.get(c.get("severity", "INFO"), "#7f8c8d"),
        "icon": "✅" if c.get("passed", False) else "❌", "detail": c.get('detail', ''),
    } for c in checks]

    return render_partial("leakage_detail.html.j2",
        safe_color="#27ae60" if leakage_audit.get('pipeline_safe') else "#e74c3c",
        safe_label="YES" if leakage_audit.get('pipeline_safe') else "NO",
        passed=leakage_audit.get('passed', 0), total_checks=leakage_audit.get('total_checks', 0),
        warnings=leakage_audit.get('warnings', 0), errors=leakage_audit.get('errors', 0), rows=rows)


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

    return render_partial("methodology.html.j2",
        pipeline_ver=pipeline_ver, formula=formula, method=method,
        n_variants=n_variants, n_traits=n_traits, n_ref=n_ref, n_pcs=n_pcs)


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

    def row_ctx(e):
        z = e.get("z_score", 0)
        pctl = round(e.get("percentile", 50), 1)
        reliable = e.get("reliable", True)
        pgs_id = e.get("pgs_id", "")
        significance = "High risk" if z>2 else ("Elevated" if z>1 else ("Low/Protective" if z<-1 else "Population average"))

        cov = pgs_coverage.get(pgs_id, {}) if pgs_coverage else {}
        n_used = cov.get("n_used", 0)
        n_total = cov.get("n_total", 0)
        snp_bar = snp_coverage_bar(n_used, n_total) if n_total > 0 else '<span style="color:#95a5a6;font-size:0.7rem">—</span>'

        bar_pct = max(5, min(95, pctl))

        return {
            "icon": "🔴" if z>2 else ("🟠" if z>1 else ("🟢" if z<-1 else "🟡")),
            "trait": e['trait'], "pgs_id": e['pgs_id'],
            "reliable_badge": reliable_badge(reliable), "snp_bar": snp_bar,
            "z_color": '#e74c3c' if z>2 else ('#f39c12' if z>1 else ('#27ae60' if z<-1 else '#2c3e50')),
            "z": f"{z:+.1f}", "pctl": pctl, "significance": significance,
            "n_snps": f"{e.get('n_snps',0):,}", "risk_bar": risk_bar(bar_pct, z),
        }

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

    # Build detailed interpretations
    detail_parts = []
    for e in high_risk + elevated:
        z = e.get('z_score', 0)
        detail_parts.append({
            "header_bg": '#fadbd8' if z>2 else '#fdebd0', "trait": e['trait'], "pgs_id": e['pgs_id'],
            "badge_class": 'high' if z>2 else 'medium',
            "badge_label": 'HIGHER RISK' if z>2 else 'ELEVATED RISK',
            "z_color": '#e74c3c' if z>2 else '#f39c12', "z": f"{z:+.1f}",
            "pctl": round(e.get('percentile',50),1), "n_snps": f"{e.get('n_snps',0):,}",
            "interpretation": interpret(e),
        })

    return render_partial("pgs_calibration.html.j2",
        total_scores=summary.get('total_scores', 0), reliable_scores=summary.get('reliable_scores', 0),
        reference_panel=methodology.get('reference_panel', '1000G')[:20],
        populations_list=', '.join(methodology.get('populations', [])),
        high_risk=[row_ctx(e) for e in high_risk], elevated=[row_ctx(e) for e in elevated],
        low_risk=[row_ctx(e) for e in low_risk], summary_html=summary_html, detail_parts=detail_parts)


def build_calibration_detail_section(calibration_report):
    """Population calibration methodology and risk breakdown."""
    methodology = calibration_report.get("methodology", {})
    thresholds = methodology.get("risk_thresholds", {})
    high_traits = calibration_report.get("high_risk_traits", [])
    medium_traits = calibration_report.get("medium_risk_traits", [])
    low_traits = calibration_report.get("low_risk_traits", [])
    populations = methodology.get("population_strata", [])

    return render_partial("calibration_detail.html.j2",
        reference_panel=methodology.get('reference_panel', 'N/A'),
        normalization=methodology.get('normalization', 'N/A').replace('_', ' ').title(),
        n_populations=len(populations), populations_list=', '.join(populations),
        traits_analyzed=calibration_report.get('traits_analyzed', 0),
        threshold_high=thresholds.get('high', '>75th'), threshold_medium=thresholds.get('medium', '25-75th'),
        threshold_low=thresholds.get('low', '<25th'),
        assigned_population=calibration_report.get('assigned_population', 'EUR'),
        calibration_note=calibration_report.get('calibration_note', '')[:300],
        n_high=len(high_traits), high_traits=high_traits,
        n_low=len(low_traits), low_traits=low_traits)


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

    def build_grouped_rows(variant_list, limit):
        """Group variants by body_system (IMPROVEMENT_PLAN.md 1.3) within a
        confidence tier, largest group first, so related findings (e.g. all
        cardiovascular hits) read together instead of in raw file order."""
        shown = variant_list[:limit]
        by_system = {}
        for v in shown:
            key = v.get("body_system") or "other"
            by_system.setdefault(key, []).append(v)
        ordered_keys = sorted(
            by_system.keys(),
            key=lambda k: (k == "other", -len(by_system[k]), BODY_SYSTEM_ORDER.index(k) if k in BODY_SYSTEM_ORDER else 99),
        )
        parts = []
        for key in ordered_keys:
            group = by_system[key]
            parts.append(
                f'<tr style="background:var(--color-bg-secondary,#f4f6f7)">'
                f'<td colspan="8" style="font-weight:700;font-size:0.8rem;padding:6px 12px">'
                f'{system_label(key, lang)} ({len(group)})</td></tr>'
            )
            parts.append("".join(build_row(v) for v in group))
        return "".join(parts)

    high_rows_html = build_grouped_rows(high_mod, 200) or '<tr><td colspan="8" style="color:var(--color-text-secondary);text-align:center;padding:1rem">✅ No high-confidence pathogenic variants found. This is normal.</td></tr>'
    low_rows_html = build_grouped_rows(low_variants, 300) or '<tr><td colspan="8" style="color:var(--color-text-secondary);text-align:center;padding:1rem">No lower-confidence variants.</td></tr>'

    return render_partial("clinvar.html.j2",
        t=t, tier_badge_high=tier_badge('high'), tier_badge_moderate=tier_badge('moderate'),
        tier_badge_low=tier_badge('low'), tier_badge_very_low=tier_badge('very_low'),
        n_high_conf=n_high_conf, n_low_conf=total - n_high_conf,
        clinvar_date=clinvar_date or 'N/A', user_vcf_total_variants=f"{meta.get('user_vcf_total_variants', 0):,}",
        n_risk=n_risk, n_high_mod=len(high_mod), high_rows_html=high_rows_html,
        n_low_variants=len(low_variants), low_rows_html=low_rows_html,
        more_not_shown=len(low_variants) - 300 if len(low_variants) > 300 else None)


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

    rows = []
    for f in findings:
        rec = f.get(f"recommendation_{lang}", f.get("recommendation_en", ""))
        if len(rec) > 150:
            rec = rec[:147] + "..."
        rows.append({
            "gene": f["gene"], "rsid": f["rsid"], "star_allele": f.get("star_allele", ""),
            "drug": f["drug"], "drug_class": f.get("drug_class", ""), "phenotype": f["phenotype"],
            "copies": f["copies"], "copies_label": "s" if f["copies"] > 1 else "",
            "recommendation": rec, "cpic_level": f["cpic_level"],
        })

    return render_partial("pharmgkb.html.j2",
        t=t, n_critical=n_critical, n_important=n_important, n_info=n_info,
        n_findings=len(findings), rows=rows)


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

def build_clinical_actionability_section(clinvar_data, pharmgkb_data, prs_entries, ui, polarity_inverted=None):
    """Cross-reference ClinVar + PharmGKB + high-risk PRS into a unified clinical summary.

    Three subsections:
    (a) High-Confidence Clinical Findings (ClinVar tier≥moderate + PharmGKB actionable)
    (b) PRS-Gene Convergence (per-trait supporting gene evidence)
    (c) Drug-Gene-PRS Intersections
    """
    lang = ui.get("_lang", "en")
    is_en = lang == "en"
    polarity_inverted = polarity_inverted or set()

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

        inverted = trait.lower() in polarity_inverted
        z_color = risk_color(z, inverted=inverted)
        if inverted:
            risk_label = ui["favorable_high"]
            label_bg, label_fg = "#d5f5e3", "#1e8449"
        else:
            risk_label = "HIGHER RISK" if is_en else "RIESGO ELEVADO"
            label_bg, label_fg = "#fadbd8", "#c0392b"

        conv_rows += (
            f'<tr><td><strong>{trait}</strong></td>'
            f'<td style="color:{z_color};font-weight:700">{z:+.2f}</td>'
            f'<td><span style="background:{label_bg};color:{label_fg};padding:2px 8px;border-radius:4px;font-size:0.75rem;font-weight:700">{risk_label}</span></td>'
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
                           ancestry=data["ancestry"],
                           polarity_inverted=data.get("_polarity_inverted", set())),
        open_by_default=True)

    # 1. Summary (always open)
    sections_html += collapsible_section("summary", f"📊 {s['summary']}",
        build_summary_cards(data["prs_result"], data["ancestry"],
                           data["integrity"], data["validation"], ui,
                           cal_lookup=data.get("_cal_lookup", {}),
                           uncert_lookup=data.get("_uncert_lookup", {}),
                           evidence_lookup=data.get("_evidence_lookup", {}),
                           portability=data.get("portability", {}),
                           polarity_inverted=data.get("_polarity_inverted", set())),
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
                       evidence_lookup=data.get("_evidence_lookup", {}),
                       polarity_inverted=data.get("_polarity_inverted", set())),
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
            ui, polarity_inverted=data.get("_polarity_inverted", set())))

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
                             evidence_lookup=data.get("_evidence_lookup"),
                             polarity_inverted=data.get("_polarity_inverted", set())))

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

    return render_document(
        lang, data, sample_id, sections_html,
        reference_coverage_banner(data["prs_result"], lang), ui,
    )


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

    data = load_report_data(args)

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
