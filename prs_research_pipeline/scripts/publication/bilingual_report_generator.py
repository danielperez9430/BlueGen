#!/usr/bin/env python3
"""
Bilingual HTML/PDF Report Generator — EN/ES.

Generates interactive HTML reports with PRS results, ancestry,
uncertainty analysis, and nutrigenetic trait interpretations.

Input:
    --prs-calibrated     CSV with columns: trait, z_score_population, percentile_population, risk_category
    --interpretations    JSON with bilingual trait interpretations
    --sample-id          Sample identifier (default: SAMPLE_001)
    --output-dir         Output directory (default: reports/)
    --lang               Language: en, es, or both (default: both)
    --ld-r2              LD threshold for methodology section (default: 0.2)

Output:
    reports/report_en.html    English HTML report
    reports/report_es.html    Spanish HTML report
    reports/report_en.pdf     English PDF (if weasyprint installed)
    reports/report_es.pdf     Spanish PDF (if weasyprint installed)

Usage:
    python bilingual_report_generator.py \\
        --prs-calibrated prs/population_calibrated_v2.csv \\
        --interpretations interpretations/interpretations.json \\
        --sample-id SAMPLE_01 --output-dir reports/ --lang both
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from html import escape as h
from pathlib import Path

import pandas as pd


# ── Constants ─────────────────────────────────────────────────────────────────

ICONS = {
    "Caffeine metabolism": "\u2615",
    "Lipid metabolism": "\u2764\uFE0F",
    "Omega-3 metabolism": "\U0001F41F",
    "Folate & methylation": "\U0001F9EC",
    "Lactose intolerance": "\U0001F95B",
    "Obesity predisposition": "\u2696\uFE0F",
    "Dopamine regulation": "\U0001F9E0",
    "Detoxification": "\U0001F6E1\uFE0F",
    "Vitamin D metabolism": "\u2600\uFE0F",
    "Glucose metabolism": "\U0001F36C",
}

COLORS = {
    "low": {"bar": "#27ae60", "bg": "#d5f5e3", "text": "#1e8449"},
    "medium": {"bar": "#f39c12", "bg": "#fdebd0", "text": "#b7950b"},
    "high": {"bar": "#e74c3c", "bg": "#fadbd8", "text": "#c0392b"},
}

POP_NAMES = {
    "en": {
        "EUR": "European",
        "AFR": "African",
        "EAS": "East Asian",
        "SAS": "South Asian",
        "AMR": "Admixed American",
    },
    "es": {
        "EUR": "Europea",
        "AFR": "Africana",
        "EAS": "Asiática Oriental",
        "SAS": "Sudasiática",
        "AMR": "Americana Mixta",
    },
}

RISK_LABELS = {
    "en": {"low": "Lower Risk", "medium": "Average Risk", "high": "Higher Risk"},
    "es": {"low": "Riesgo Bajo", "medium": "Riesgo Promedio", "high": "Riesgo Elevado"},
}

CSS = """*{box-sizing:border-box;margin:0;padding:0}body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;color:#2c3e50;background:#f8f9fa;line-height:1.6}.report-header{background:linear-gradient(135deg,#1a1a2e,#16213e,#0f3460);color:#fff;padding:3rem 2rem;text-align:center}.container{max-width:960px;margin:0 auto;padding:2rem 1.5rem}.risk-category-badge{display:inline-block;padding:.15rem .6rem;border-radius:4px;font-size:.75rem;font-weight:700;text-transform:uppercase}.risk-category-badge.low{background:#d5f5e3;color:#1e8449}.risk-category-badge.medium{background:#fdebd0;color:#b7950b}.risk-category-badge.high{background:#fadbd8;color:#c0392b}table{width:100%;border-collapse:collapse}td{padding:.65rem .75rem;border-top:1px solid #dee2e6;font-size:.9rem}tr:hover{background:#f8f9fa}.disclaimer-box{background:#fef9e7;border:1px solid #f9e79f;border-radius:8px;padding:1.25rem;margin-top:2rem;font-size:.8rem;color:#7d6608}.report-footer{text-align:center;padding:2rem;font-size:.75rem;color:#7f8c8d;border-top:1px solid #dee2e6;margin-top:2rem}@media print{body{background:#fff}}"""


def _t(en: str, es: str, is_es: bool) -> str:
    return es if is_es else en


def _load_json(path: Path) -> dict:
    if path.exists():
        with open(path) as fh:
            return json.load(fh)
    return {}


def _load_csv(path: Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path)
    return None


def _find_ancestry() -> dict:
    for p in ["ancestry/classification_report.json", "pca/ancestry_inference.json"]:
        data = _load_json(Path(p))
        if data:
            return data
    return {}


def generate_report(
    prs_calibrated: str,
    interpretations: str,
    sample_id: str = "SAMPLE_001",
    output_dir: str = "reports/",
    lang: str = "both",
    ld_r2: float = 0.2,
) -> list[str]:
    """Generate bilingual HTML reports. Returns list of generated file paths."""
    prs = pd.read_csv(prs_calibrated)
    interps = _load_json(Path(interpretations))
    ancestry = _find_ancestry()
    uncertainty = _load_csv(Path("prs/prs_uncertainty.csv"))
    consistency = _load_json(Path("prs/consistency_check_report.json"))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    generated = []

    langs = ["en", "es"] if lang == "both" else [lang]

    for lg in langs:
        is_es = lg == "es"
        interp = interps.get(lg, interps.get("en", {}))
        meta = interp.get("metadata", {})
        pop = meta.get("assigned_population", "EUR")
        pop_name = POP_NAMES[lg].get(pop, pop)
        hi = meta.get("high_risk_count", 0)
        md = meta.get("medium_risk_count", 0)
        lo = meta.get("low_risk_count", 0)
        date_str = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

        # PRS table rows
        prs_rows = ""
        for _, row in prs.iterrows():
            trait = row["trait"]
            risk = row.get("risk_category", "medium")
            z = float(row.get("z_score_population", 0))
            pctl = float(row.get("percentile_population", 50))
            icon = ICONS.get(trait, "\U0001F9EC")
            c = COLORS.get(risk, COLORS["medium"])
            bar = int(min(max((z + 3) / 6 * 100, 0), 100))
            prs_rows += (
                f'<tr><td>{icon} {h(trait)}</td>'
                f'<td><strong>{z:+.2f}</strong></td>'
                f'<td>{pctl:.1f}%</td>'
                f'<td><span class="risk-category-badge {risk}">'
                f'{RISK_LABELS[lg].get(risk, risk).upper()}</span></td>'
                f'<td><div style="display:flex;align-items:center">'
                f'<div style="flex:1;height:10px;background:#e9ecef;border-radius:5px;overflow:hidden">'
                f'<div style="width:{bar}%;height:100%;background:{c["bar"]};border-radius:5px">'
                f"</div></div></div></td></tr>"
            )

        # Trait detail sections
        ti_map = interp.get("trait_interpretations", {})
        trait_secs = ""
        for _, row in prs.iterrows():
            trait = row["trait"]
            risk = row.get("risk_category", "medium")
            z = float(row.get("z_score_population", 0))
            pctl = float(row.get("percentile_population", 50))
            icon = ICONS.get(trait, "\U0001F9EC")
            c = COLORS.get(risk, COLORS["medium"])
            bar = int(min(max((z + 3) / 6 * 100, 0), 100))
            ti = ti_map.get(trait, {})
            dietary = "".join(
                f"<li>{h(item)}</li>" for item in ti.get("dietary_context", [])
            )
            dietary_html = (
                f'<div style="background:#eaf2f8;border-radius:8px;padding:1rem;margin-top:1rem">'
                f"<h4>{_t('Dietary & Lifestyle Context', 'Contexto Dietético', is_es)}</h4>"
                f"<ul>{dietary}</ul></div>"
                if dietary
                else ""
            )
            trait_secs += (
                f'<div style="background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08);'
                f'margin-bottom:1.5rem;overflow:hidden">'
                f'<div style="padding:1rem 1.25rem;display:flex;justify-content:space-between;'
                f'align-items:center;background:{c["bg"]}">'
                f'<span style="font-size:1.1rem;font-weight:600">{icon} {h(trait)}</span>'
                f'<span><span class="risk-category-badge {risk}">'
                f'{RISK_LABELS[lg].get(risk, risk).upper()}</span></span></div>'
                f'<div style="padding:1.25rem;border-top:1px solid #dee2e6">'
                f'<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(110px,1fr));'
                f'gap:1rem;margin-bottom:1rem">'
                f'<div style="text-align:center">'
                f'<div style="font-size:1.3rem;font-weight:700;color:{c["text"]}">{z:+.2f}</div>'
                f'<div style="font-size:.75rem;color:#7f8c8d">Z-Score</div></div>'
                f'<div style="text-align:center">'
                f'<div style="font-size:1.3rem;font-weight:700">{pctl:.1f}%</div>'
                f'<div style="font-size:.75rem;color:#7f8c8d">'
                f'{_t("Percentile", "Percentil", is_es)}</div></div></div>'
                f'<div style="display:flex;align-items:center;gap:.5rem;margin:1rem 0">'
                f'<span style="font-size:.7rem;color:#27ae60">'
                f'{_t("Low", "Bajo", is_es)}</span>'
                f'<div style="flex:1;height:10px;background:#e9ecef;border-radius:5px;overflow:hidden">'
                f'<div style="width:{bar}%;height:100%;background:{c["bar"]};border-radius:5px">'
                f'</div></div>'
                f'<span style="font-size:.7rem;color:#e74c3c">'
                f'{_t("High", "Alto", is_es)}</span></div>'
                f'<div style="background:#f8f9fa;border-radius:8px;padding:1rem;margin:1rem 0;'
                f'font-size:.9rem;line-height:1.7">'
                f"<strong>\U0001F9EC {h(ti.get('gene', ''))}</strong>"
                f'<p style="color:#7f8c8d;font-size:.8rem;margin-top:.25rem">'
                f"{h(ti.get('pathway', ''))}</p>"
                f'<p style="margin-top:.75rem;white-space:pre-line">'
                f"{h(ti.get('description', ''))}</p>"
                f'<p style="margin-top:.75rem;white-space:pre-line;font-weight:500">'
                f"{h(ti.get('narrative', ''))}</p></div>{dietary_html}</div></div>"
            )

        # Uncertainty section
        unc_html = ""
        if uncertainty is not None and len(uncertainty) > 0:
            u_rows = "".join(
                f'<tr><td>{h(r["trait"])}</td>'
                f'<td>{float(r["prs"]):.3f} ± {float(r["prs_se"]):.3f}</td>'
                f'<td>[{float(r["ci_95_lower"]):.3f}, {float(r["ci_95_upper"]):.3f}]</td>'
                f'<td>{float(r["uncertainty_score"]):.3f}</td></tr>'
                for _, r in uncertainty.iterrows()
            )
            unc_html = (
                f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
                f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
                f'📊 {_t("Uncertainty Overview", "Resumen de Incertidumbre", is_es)}</h2>'
                f'<table style="width:100%;border-collapse:collapse;background:#fff;border-radius:8px;'
                f'overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);margin:1rem 0">'
                f'<thead><tr style="background:#f1f3f5">'
                f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
                f'text-transform:uppercase">'
                f'{_t("Trait", "Rasgo", is_es)}</th>'
                f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
                f'text-transform:uppercase">PRS ± SE</th>'
                f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
                f'text-transform:uppercase">95% CI</th>'
                f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
                f'text-transform:uppercase">'
                f'{_t("Uncert.", "Incert.", is_es)}</th>'
                f"</tr></thead><tbody>{u_rows}</tbody></table></section>"
            )

        # Consistency section
        cons_html = ""
        if consistency:
            passed = consistency.get("passed", True)
            sc = "#27ae60" if passed else "#e74c3c"
            st = (
                ("\u2705 " + _t("PASSED", "COMPATIBLE", is_es))
                if passed
                else ("\u274C " + _t("FAILED", "INCOMPATIBLE", is_es))
            )
            cons_html = (
                f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
                f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
                f'🔬 {_t("GWAS-Ancestry Compatibility", "Compatibilidad GWAS-Ancestral", is_es)}</h2>'
                f'<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));'
                f'gap:1rem;margin:1.5rem 0">'
                f'<div style="background:#fff;border-radius:8px;padding:1.25rem;text-align:center;'
                f'box-shadow:0 1px 3px rgba(0,0,0,.08)">'
                f'<div style="font-size:2rem;font-weight:800;color:{sc}">{st}</div>'
                f'<div style="font-size:.85rem;color:#7f8c8d">'
                f'{_t("Validation Status", "Estado de Validación", is_es)}</div></div></div></section>'
            )

        # Full HTML
        html = (
            f'<!DOCTYPE html><html lang="{lg}">'
            f'<head><meta charset="UTF-8"><title>'
            f'{_t("PRS Research Report", "Informe PRS", is_es)} — {h(sample_id)}</title>'
            f"<style>{CSS}</style></head><body>"
            f'<header class="report-header">'
            f'<h1 style="font-size:2.2rem">'
            f'{_t("PRS Research Report", "Informe de Investigación PRS", is_es)}</h1>'
            f'<div style="font-size:1.1rem;opacity:.9">'
            f'{_t("Population-Calibrated PRS Analysis", "Análisis de PRS Calibrado por Población", is_es)}</div>'
            f'<div style="margin-top:1rem;font-size:.85rem;opacity:.7">'
            f'{_t("Sample", "Muestra", is_es)}: {h(sample_id)} | '
            f'{_t("Date", "Fecha", is_es)}: {date_str} | Pipeline v6.0.0 | GRCh37/hg19</div></header>'
            f'<div class="container">'
            # Genetic Overview
            f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
            f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
            f'📊 {_t("Genetic Overview", "Resumen Genético", is_es)}</h2>'
            f'<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));'
            f'gap:1rem;margin-bottom:2rem">'
            f'<div style="background:#fff;border-radius:8px;padding:1.25rem;text-align:center;'
            f'box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #e74c3c">'
            f'<div style="font-size:2.2rem;font-weight:800;color:#e74c3c">{hi}</div>'
            f'<div style="font-size:.85rem;color:#7f8c8d">'
            f'{_t("Higher Risk", "Riesgo Elevado", is_es)}</div></div>'
            f'<div style="background:#fff;border-radius:8px;padding:1.25rem;text-align:center;'
            f'box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #f39c12">'
            f'<div style="font-size:2.2rem;font-weight:800;color:#f39c12">{md}</div>'
            f'<div style="font-size:.85rem;color:#7f8c8d">'
            f'{_t("Average Risk", "Riesgo Promedio", is_es)}</div></div>'
            f'<div style="background:#fff;border-radius:8px;padding:1.25rem;text-align:center;'
            f'box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #27ae60">'
            f'<div style="font-size:2.2rem;font-weight:800;color:#27ae60">{lo}</div>'
            f'<div style="font-size:.85rem;color:#7f8c8d">'
            f'{_t("Lower Risk", "Riesgo Bajo", is_es)}</div></div></div>'
            f'<p style="white-space:pre-line;background:#fff;padding:1.25rem;border-radius:8px;'
            f'box-shadow:0 1px 3px rgba(0,0,0,.08)">'
            f'{h(interp.get("global_summary", ""))}</p></section>'
            # Ancestry
            f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
            f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
            f'\U0001F30D {_t("Ancestry", "Ascendencia", is_es)}</h2>'
            f'<p style="font-size:1.1rem"><strong>{pop_name} ({pop})</strong></p></section>'
            # PRS Table
            f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
            f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
            f'📈 {_t("Population-Calibrated PRS", "PRS Calibrado por Población", is_es)}</h2>'
            f'<table><thead><tr style="background:#f1f3f5">'
            f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
            f'text-transform:uppercase">{_t("Trait", "Rasgo", is_es)}</th>'
            f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
            f'text-transform:uppercase">Z-Score</th>'
            f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
            f'text-transform:uppercase">{_t("Percentile", "Percentil", is_es)}</th>'
            f'<th style="padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;'
            f'text-transform:uppercase">{_t("Risk", "Riesgo", is_es)}</th>'
            f'</tr></thead><tbody>{prs_rows}</tbody></table></section>'
            # Trait Analysis
            f'<h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
            f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
            f'🔬 {_t("Nutrigenetic Trait Analysis", "Análisis Nutrigenético", is_es)}</h2>'
            f"{trait_secs}{unc_html}{cons_html}"
            # Methodology
            f'<section><h2 style="font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;'
            f'padding-bottom:.5rem;border-bottom:2px solid #dee2e6">'
            f'📐 {_t("Methodology", "Metodología", is_es)}</h2>'
            f'<p><strong>{_t("PRS Model:", "Modelo PRS:", is_es)}</strong> '
            f"PRS = Σ (βₜ × Gᵢⱼ) "
            f'{_t("with", "con", is_es)} LD pruning (r² &lt; {ld_r2}), '
            f"true PCA projection onto 1000 Genomes, empirical population calibration, "
            f"3-layer uncertainty propagation, Phase 6 corrections.</p></section>"
            # Disclaimer
            f'<div class="disclaimer-box"><h3 style="color:#9a7d0a;margin-bottom:.5rem">'
            f'\u26A0\uFE0F {_t("Important Disclaimer", "Aviso Importante", is_es)}</h3>'
            f'<p style="white-space:pre-line">'
            f'{h(interp.get("disclaimer", ""))}</p></div>'
            f"</div>"
            # Footer
            f'<footer class="report-footer"><p>'
            f'{_t("PRS Research Platform", "Plataforma PRS", is_es)} v6.0.0 — '
            f'{_t("Generated", "Generado", is_es)} {date_str}</p>'
            f'<p>{_t("Sample", "Muestra", is_es)}: {h(sample_id)} | GRCh37/hg19</p></footer>'
            f"</body></html>"
        )

        html_path = out / f"report_{lg}.html"
        with open(html_path, "w", encoding="utf-8") as fh:
            fh.write(html)
        generated.append(str(html_path))

        try:
            from weasyprint import HTML

            pdf_path = out / f"report_{lg}.pdf"
            HTML(string=html).write_pdf(str(pdf_path))
            generated.append(str(pdf_path))
            print(f"\u2705 [{lg}] {html_path} + {pdf_path}")
        except Exception as e:
            print(f"\u2705 [{lg}] {html_path}  \u26A0\uFE0F PDF: {e}")

    print("\u2705 Bilingual reports generated (v6.0.0)")
    return generated


def main():
    parser = argparse.ArgumentParser(description="Bilingual HTML/PDF Report Generator")
    parser.add_argument(
        "--prs-calibrated", required=True, help="Path to calibrated PRS CSV"
    )
    parser.add_argument(
        "--interpretations", required=True, help="Path to interpretations JSON"
    )
    parser.add_argument("--sample-id", default="SAMPLE_001", help="Sample identifier")
    parser.add_argument("--output-dir", default="reports/", help="Output directory")
    parser.add_argument(
        "--lang", default="both", choices=["en", "es", "both"], help="Report language"
    )
    parser.add_argument("--ld-r2", type=float, default=0.2, help="LD r2 threshold")
    args = parser.parse_args()

    generate_report(
        prs_calibrated=args.prs_calibrated,
        interpretations=args.interpretations,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
        lang=args.lang,
        ld_r2=args.ld_r2,
    )


if __name__ == "__main__":
    main()
