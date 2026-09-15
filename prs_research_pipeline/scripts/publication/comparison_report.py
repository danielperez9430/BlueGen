#!/usr/bin/env python3
"""
Multi-sample comparison report (RELEASE_PLAN 3.0.4).

`prs.py run --vcf a.vcf.gz,b.vcf.gz` (or a multi-sample VCF) already scores
every sample: Stage F/H write one row per (individual, trait) and the PGS
calibration one per (individual, score). The comprehensive report is
single-sample, so this builds reports/comparison_report_{en,es}.html with
the samples side by side:

  - curated PRS: population z-score / percentile / risk per trait,
  - PGS Catalog: z-score / percentile per score,
  - ancestry: the assigned super-population per sample when the classifier
    wrote per-sample results (otherwise the run-level one, marked as such).

No relatedness or inheritance inference: it is a table, not a family
analysis. ClinVar and pharmacogenomics stay in the per-sample report.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report.render import _env, _read_static  # noqa: E402
from utils.constants import PIPELINE_VERSION  # noqa: E402

UI = {
    "en": {"title": "BlueGen — Sample Comparison", "subtitle": "Side-by-side polygenic scores for a multi-sample run",
           "prs": "Curated PRS (population-calibrated)", "pgs": "PGS Catalog (calibrated)",
           "ancestry": "Ancestry", "trait": "Trait", "score": "Score", "no_prs": "No calibrated PRS found.",
           "no_pgs": "No calibrated PGS Catalog scores found.", "legend": "Cells show z-score (percentile). "
           "Colour: red ≥ +1σ, green ≤ −1σ, grey in between; ‡ = low-confidence reference distribution.",
           "note": "Research use only. ClinVar and pharmacogenomics are reported per sample in the comprehensive report. "
                   "No relatedness is inferred here.", "samples": "samples", "run_level": "run-level (not per sample)"},
    "es": {"title": "BlueGen — Comparación de muestras", "subtitle": "Scores poligénicos lado a lado de una ejecución multi-muestra",
           "prs": "PRS curado (calibrado por población)", "pgs": "PGS Catalog (calibrado)",
           "ancestry": "Ancestría", "trait": "Rasgo", "score": "Score", "no_prs": "No hay PRS calibrados.",
           "no_pgs": "No hay scores del PGS Catalog calibrados.", "legend": "Cada celda muestra z (percentil). "
           "Color: rojo ≥ +1σ, verde ≤ −1σ, gris entre medias; ‡ = distribución de referencia de baja confianza.",
           "note": "Solo uso en investigación. ClinVar y farmacogenómica se informan por muestra en el informe completo. "
                   "Aquí no se infiere parentesco.", "samples": "muestras", "run_level": "de la ejecución (no por muestra)"},
}


def _cell(z, pctl, low_conf=False):
    if z is None or pd.isna(z):
        return {"text": "—", "cls": "na"}
    cls = "hi" if z >= 1 else ("lo" if z <= -1 else "mid")
    mark = " ‡" if low_conf else ""
    return {"text": f"{z:+.2f} ({pctl:.0f}%){mark}", "cls": cls}


def load_tables(pipeline_dir: Path):
    prs_path = pipeline_dir / "prs" / "population_calibrated_v2.csv"
    if not prs_path.exists():
        prs_path = pipeline_dir / "prs" / "population_calibrated.csv"
    prs = pd.read_csv(prs_path) if prs_path.exists() else pd.DataFrame()
    pgs_path = pipeline_dir / "prs" / "pgs_scores" / "pgs_calibrated.csv"
    pgs = pd.read_csv(pgs_path) if pgs_path.exists() else pd.DataFrame()
    anc_path = pipeline_dir / "science" / "ANCESTRY_MODEL.json"
    ancestry = json.loads(anc_path.read_text()) if anc_path.exists() else {}
    return prs, pgs, ancestry


def build_context(prs: pd.DataFrame, pgs: pd.DataFrame, ancestry: dict, lang: str) -> dict:
    ui = UI[lang]
    samples = []
    for df in (prs, pgs):
        if not df.empty and "individual_id" in df.columns:
            for s in df["individual_id"].astype(str).unique():
                if s not in samples:
                    samples.append(s)

    prs_rows = []
    if not prs.empty and {"individual_id", "trait"} <= set(prs.columns):
        zc = "z_score_population" if "z_score_population" in prs.columns else "population_zscore"
        pc = "percentile_population" if "percentile_population" in prs.columns else "population_percentile"
        lc = "low_confidence" if "low_confidence" in prs.columns else None
        for trait in sorted(prs["trait"].astype(str).unique()):
            sub = prs[prs["trait"].astype(str) == trait]
            cells = []
            for s in samples:
                r = sub[sub["individual_id"].astype(str) == s]
                if r.empty:
                    cells.append(_cell(None, None))
                else:
                    r = r.iloc[0]
                    cells.append(_cell(float(r[zc]), float(r[pc]), bool(r[lc]) if lc else False))
            spread = [float(sub[sub["individual_id"].astype(str) == s][zc].iloc[0]) for s in samples
                      if not sub[sub["individual_id"].astype(str) == s].empty]
            prs_rows.append({"name": trait, "cells": cells,
                             "spread": (max(spread) - min(spread)) if len(spread) > 1 else 0.0})
        prs_rows.sort(key=lambda r: -r["spread"])

    pgs_rows = []
    if not pgs.empty and {"individual_id", "pgs_id"} <= set(pgs.columns):
        for pgs_id in sorted(pgs["pgs_id"].astype(str).unique()):
            sub = pgs[pgs["pgs_id"].astype(str) == pgs_id]
            trait = str(sub["trait"].iloc[0]) if "trait" in sub.columns else ""
            cells = []
            for s in samples:
                r = sub[sub["individual_id"].astype(str) == s]
                if r.empty:
                    cells.append(_cell(None, None))
                else:
                    r = r.iloc[0]
                    rel = ("reliable" in r.index) and not bool(r["reliable"])
                    cells.append(_cell(float(r["z_score"]), float(r["percentile"]), rel))
            spread = [float(sub[sub["individual_id"].astype(str) == s]["z_score"].iloc[0]) for s in samples
                      if not sub[sub["individual_id"].astype(str) == s].empty]
            pgs_rows.append({"name": f"{trait} ({pgs_id})" if trait else pgs_id, "cells": cells,
                             "spread": (max(spread) - min(spread)) if len(spread) > 1 else 0.0})
        pgs_rows.sort(key=lambda r: -r["spread"])

    per_sample = ancestry.get("per_sample") if isinstance(ancestry, dict) else None
    if isinstance(per_sample, dict) and per_sample:
        anc_cells = [str(per_sample.get(s, {}).get("assigned_population", "—")) for s in samples]
        anc_note = ""
    else:
        pop = ancestry.get("assigned_population", "—") if isinstance(ancestry, dict) else "—"
        anc_cells = [str(pop)] * len(samples)
        anc_note = ui["run_level"]

    return {
        "lang": lang, "ui": ui, "samples": samples, "n_samples": len(samples),
        "prs_rows": prs_rows, "pgs_rows": pgs_rows, "anc_cells": anc_cells, "anc_note": anc_note,
        "now": datetime.now().strftime("%Y-%m-%d %H:%M UTC"), "pipeline_version": PIPELINE_VERSION,
        "css": _read_static("report.css"),
    }


def render(prs, pgs, ancestry, lang: str) -> str:
    return _env.get_template("comparison.html.j2").render(**build_context(prs, pgs, ancestry, lang))


def main(argv=None):
    p = argparse.ArgumentParser(description="Multi-sample comparison report")
    p.add_argument("--pipeline-dir", default=".", help="prs_research_pipeline directory (default: cwd)")
    p.add_argument("--output-dir", default="reports/")
    p.add_argument("--lang", default="both", choices=["en", "es", "both"])
    args = p.parse_args(argv)

    pipeline_dir = Path(args.pipeline_dir)
    prs, pgs, ancestry = load_tables(pipeline_dir)
    n = len(set(prs.get("individual_id", pd.Series(dtype=str)).astype(str))
            | set(pgs.get("individual_id", pd.Series(dtype=str)).astype(str)))
    if n < 2:
        print(f"  Comparison report: {n} sample(s) — nothing to compare")
        return 0
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for lang in (["en", "es"] if args.lang == "both" else [args.lang]):
        path = out_dir / f"comparison_report_{lang}.html"
        path.write_text(render(prs, pgs, ancestry, lang), encoding="utf-8")
        print(f"  ✓ {path} ({n} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
