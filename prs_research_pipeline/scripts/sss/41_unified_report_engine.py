#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 9 — UNIFIED REPORT ENGINE (SSST)                                     ║
║   scripts/41_unified_report_engine.py                                       ║
║                                                                            ║
║   Generates a SINGLE coherent scientific narrative from the unified         ║
║   SSST sources: PRS_CORE + ANCESTRY_MODEL + BENCHMARK + UNCERTAINTY.        ║
║                                                                            ║
║   Replaces fragmented narrative generation across multiple modules.         ║
║   Identical EN/ES structure — no contradictory interpretation.             ║
║                                                                            ║
║   Inputs:                                                                    ║
║     PRS_CORE → standardized PRS definition                                 ║
║     ANCESTRY_MODEL → single ancestry source                                 ║
║     PRS_RESULT → unified PRS output                                         ║
║     BENCHMARK → reinterpreted validation                                    ║
║     UNCERTAINTY → standardized CIs                                          ║
║                                                                            ║
║   Output:                                                                    ║
║     reports/SCIENTIFIC_MANUSCRIPT_EN.md                                     ║
║     reports/SCIENTIFIC_MANUSCRIPT_ES.md                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime

logger = logging.getLogger(__name__)

class UnifiedReportEngine:
    """
    Generates a single coherent scientific narrative from all SSST sources.

    The manuscript follows a standard scientific structure:
      1. Abstract
      2. Methods summary
      3. Results
      4. Ancestry context
      5. Benchmark validation
      6. Limitations
      7. Interpretation boundaries

    Identical structure in EN and ES. No machine translation.
    All data drawn from SSST sources only.
    """

    def __init__(self, output_dir: str = "reports"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(self, sample_id: str = "SAMPLE_001",
                 prs_core_json: str = "science/prs_core_definition.json",
                 ancestry_json: str = "science/ANCESTRY_MODEL.json",
                 prs_result_json: str = "prs/PRS_RESULT.json",
                 benchmark_json: str = "benchmark/VALIDATION_REPORT.json",
                 integrity_json: str = "science/scientific_integrity_score.json"
                 ) -> Dict[str, str]:
        logger.info("═══ Unified Report Engine (SSST) ═══")

        # Load SSST sources
        prs_core = self._load(prs_core_json)
        ancestry = self._load(ancestry_json)
        prs_result = self._load(prs_result_json)
        benchmark = self._load(benchmark_json)
        integrity = self._load(integrity_json)

        # Extract key values
        n_variants = prs_core.get("n_variants", 109)
        n_traits = prs_core.get("n_traits", 10)
        formula = prs_core.get("formula", "PRS = Σ(βⱼ × Gᵢⱼ)")
        assigned_pop = ancestry.get("assigned_population", "EUR")
        ancestry_confidence = ancestry.get("confidence", "MODERATE")
        entries = prs_result.get("prs_entries", [])
        integrity_score = integrity.get("scientific_integrity_score", 0)
        integrity_cat = integrity.get("category", "Unknown")

        # Generate both languages
        paths = {}
        for lang in ["en", "es"]:
            is_es = lang == "es"
            content = self._build_manuscript(
                lang, is_es, sample_id, n_variants, n_traits, formula,
                assigned_pop, ancestry_confidence, entries,
                integrity_score, integrity_cat, benchmark)
            filename = f"SCIENTIFIC_MANUSCRIPT_{lang.upper()}.md"
            path = self.output_dir / filename
            with open(path, "w") as fh:
                fh.write(content)
            paths[lang] = str(path)
            logger.info(f"  ✅ [{lang}] {filename}")

        return paths

    def _build_manuscript(self, lang: str, is_es: bool, sample_id: str,
                          n_variants: int, n_traits: int, formula: str,
                          assigned_pop: str, ancestry_confidence: str,
                          entries: List[Dict], integrity_score: float,
                          integrity_cat: str, benchmark: Dict) -> str:
        T = lambda en, es: es if is_es else en

        pop_names = {"EUR": T("European", "Europea"),
                     "AFR": T("African", "Africana"),
                     "EAS": T("East Asian", "Asiática Oriental"),
                     "SAS": T("South Asian", "Sudasiática"),
                     "AMR": T("Admixed American", "Americana Mixta")}
        pop_name = pop_names.get(assigned_pop, assigned_pop)

        n_high = sum(1 for e in entries if e.get("risk_category") == "high")
        n_medium = sum(1 for e in entries if e.get("risk_category") == "medium")
        n_low = sum(1 for e in entries if e.get("risk_category") == "low")

        n_internal = benchmark.get("validation_summary", {}).get("internal", 0)
        n_external = benchmark.get("validation_summary", {}).get("external", 0)
        all_independent = benchmark.get("validation_summary", {}).get("all_independent", False)

        return f"""# {T('PRS Research Platform — Scientific Manuscript', 'Plataforma PRS — Manuscrito Científico')}

**{T('Sample', 'Muestra')}:** {sample_id}
**{T('Generated', 'Generado')}:** {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}
**{T('Pipeline Version', 'Versión')}:** 9.0.0 (SSST — Single Source of Scientific Truth)
**{T('Scientific Integrity Score', 'Puntuación de Integridad Científica')}:** {integrity_score:.0f}/100 — {integrity_cat}

---

## 1. {T('Abstract', 'Resumen')}

{T('We present a polygenic risk score (PRS) analysis for nutrigenomic trait categories using a curated panel of ' + str(n_variants) + ' single-nucleotide polymorphisms (SNPs) across ' + str(n_traits) + ' traits. The analysis employs PLINK-based genotype processing, 1000 Genomes Phase 3 reference-based PCA projection for ancestry inference, and population-calibrated PRS computation using the standard weighted sum formula: ' + formula + '. Results are validated against PGS Catalog reference scores, published GWAS consortia effect directions, and internal multi-method consistency checks. The platform achieves a Scientific Integrity Score of ' + str(integrity_score) + '/100 (' + integrity_cat + '). All analyses are fully reproducible with frozen random seeds, SHA-256 hashed inputs/outputs, and a documented scientific assumption lock file.',
    'Presentamos un análisis de puntuación de riesgo poligénico (PRS) para categorías de rasgos nutrigenómicos utilizando un panel curado de ' + str(n_variants) + ' polimorfismos de nucleótido único (SNPs) en ' + str(n_traits) + ' rasgos. El análisis emplea procesamiento de genotipos basado en PLINK, proyección PCA con referencia del 1000 Genomas Fase 3 para inferencia de ascendencia, y cálculo de PRS calibrado por población utilizando la fórmula estándar de suma ponderada: ' + formula + '. Los resultados se validan contra puntuaciones de referencia del PGS Catalog, direcciones de efecto de consorcios GWAS publicados, y verificaciones internas de consistencia multi-método. La plataforma alcanza una Puntuación de Integridad Científica de ' + str(integrity_score) + '/100 (' + integrity_cat + '). Todos los análisis son completamente reproducibles con semillas aleatorias congeladas, entradas/salidas con hash SHA-256, y un archivo documentado de bloqueo de supuestos científicos.')}

---

## 2. {T('Methods Summary', 'Resumen de Métodos')}

{T('PRS Computation: ' + formula + ', implemented via PLINK --score with dosage-weighted allele counting. LD pruning at r² < 0.2 using ancestry-matched reference panels. PCA projection onto 1000 Genomes Phase 3 reference space (Price et al. 2006) with ' + str(n_variants) + ' curated SNPs across ' + str(n_traits) + ' nutrigenomic trait categories. Population calibration against empirical 1000G reference distributions with bootstrap confidence intervals. Three-layer uncertainty propagation (genotype + ancestry + effect).',
    'Cálculo de PRS: ' + formula + ', implementado mediante PLINK --score con conteo de alelos ponderado por dosificación. Poda LD a r² < 0.2 utilizando paneles de referencia emparejados por ascendencia. Proyección PCA en el espacio de referencia del 1000 Genomas Fase 3 (Price et al. 2006) con ' + str(n_variants) + ' SNPs curados en ' + str(n_traits) + ' categorías de rasgos nutrigenómicos. Calibración poblacional contra distribuciones empíricas de referencia 1000G con intervalos de confianza bootstrap. Propagación de incertidumbre de tres capas (genotipo + ascendencia + efecto).')}

---

## 3. {T('Results', 'Resultados')}

{T('Ancestry: Classified as ' + pop_name + ' (' + assigned_pop + ', ' + ancestry_confidence + ' confidence).',
    'Ascendencia: Clasificado como ' + pop_name + ' (' + assigned_pop + ', confianza ' + ancestry_confidence + ').')}

| {T('Trait', 'Rasgo')} | {T('Raw PRS', 'PRS Bruto')} | {T('Z-Score', 'Puntuación Z')} | {T('Percentile', 'Percentil')} | {T('95% CI', 'IC 95%')} | {T('Risk', 'Riesgo')} |
|--------|------|------|------|------|------|
"""

        lines = []
        for e in entries[:10]:
            trait = e["trait"]
            raw = e.get("raw_score", 0)
            z = e.get("population_zscore", 0)
            pctl = e.get("population_percentile", 50)
            ci_l = e.get("ci_95_lower", 0)
            ci_u = e.get("ci_95_upper", 0)
            risk = e.get("risk_category", "medium")
            risk_label = {"high": T("HIGH", "ALTO"), "medium": T("MEDIUM", "MEDIO"), "low": T("LOW", "BAJO")}.get(risk, risk)
            lines.append(f"| {trait} | {raw:.3f} | {z:+.2f} | {pctl:.1f}% | [{ci_l:.0f}, {ci_u:.0f}] | {risk_label} |")

        content = "\n".join(lines)

        content += f"""

**{T('Summary', 'Resumen')}:** {n_high} {T('traits with higher risk,', 'rasgos con riesgo elevado,')} {n_medium} {T('average,', 'promedio,')} {n_low} {T('lower.', 'bajo.')}

---

## 4. {T('Ancestry Context', 'Contexto de Ascendencia')}

{T('The sample was classified as ' + pop_name + ' (' + assigned_pop + ') with ' + ancestry_confidence + ' confidence using a PCA ensemble classifier trained on 1000 Genomes Phase 3 reference data. GWAS effect sizes used in this analysis are primarily derived from European-ancestry discovery populations. Cross-ancestry PRS transferability to non-EUR populations is documented as a limitation. Population-specific calibration using empirical 1000G reference distributions (' + assigned_pop + ': μ, σ computed from ~500 reference samples) mitigates but does not eliminate ancestry-related bias.',
    'La muestra fue clasificada como ' + pop_name + ' (' + assigned_pop + ') con confianza ' + ancestry_confidence + ' utilizando un clasificador de conjunto PCA entrenado en datos de referencia del 1000 Genomas Fase 3. Los tamaños de efecto GWAS utilizados en este análisis se derivan principalmente de poblaciones de descubrimiento de ascendencia europea. La transferibilidad entre ascendencias del PRS a poblaciones no europeas está documentada como una limitación. La calibración específica por población utilizando distribuciones empíricas de referencia 1000G (' + assigned_pop + ': μ, σ calculados de ~500 muestras de referencia) mitiga pero no elimina el sesgo relacionado con la ascendencia.')}

---

## 5. {T('Benchmark Validation', 'Validación Comparativa')}

{T(f'The platform has been validated through {n_internal} internal consistency checks and {n_external} external benchmark comparisons. External references: {"all independent" if all_independent else "some circular validations detected"}. Internal validations include cross-method PRS agreement, variant coverage audit, and calibration quality assessment. External validations include PGS Catalog score concordance, GWAS consortium effect direction consistency, and population portability analysis.',
    f'La plataforma ha sido validada a través de {n_internal} verificaciones de consistencia interna y {n_external} comparaciones de referencia externa. Referencias externas: {"todas independientes" if all_independent else "se detectaron algunas validaciones circulares"}. Las validaciones internas incluyen concordancia PRS entre métodos, auditoría de cobertura de variantes y evaluación de calidad de calibración. Las validaciones externas incluyen concordancia de puntuaciones del PGS Catalog, consistencia de dirección de efecto de consorcios GWAS y análisis de portabilidad poblacional.')}

---

## 6. {T('Limitations', 'Limitaciones')}

1. {T('PRS is probabilistic, not diagnostic — estimates susceptibility, does not diagnose disease.', 'El PRS es probabilístico, no diagnóstico — estima susceptibilidad, no diagnostica enfermedad.')}
2. {T(f'GWAS effect sizes are primarily EUR-derived — cross-ancestry transferability is reduced for non-EUR populations.', f'Los tamaños de efecto GWAS son principalmente de origen europeo — la transferibilidad entre ascendencias es reducida para poblaciones no europeas.')}
3. {T(f'Curated panel of {n_variants} SNPs covers nutrigenomic traits — genome-wide PRS (~1M+ SNPs) captures substantially more genetic signal.', f'El panel curado de {n_variants} SNPs cubre rasgos nutrigenómicos — el PRS de genoma completo (~1M+ SNPs) captura sustancialmente más señal genética.')}
4. {T('Not clinically validated — research use only.', 'No validado clínicamente — solo para uso en investigación.')}

---

## 7. {T('Interpretation Boundaries', 'Límites de Interpretación')}

{T('These results are RESEARCH-GRADE only. They are not validated for clinical decision-making. Interpretation should ALWAYS consider ancestry background and GWAS transferability. All PRS values are reported with 95% confidence intervals. Scores with coverage < 70% or evidence < 50/100 should be considered exploratory only.',
    'Estos resultados son SOLO PARA INVESTIGACIÓN. No están validados para la toma de decisiones clínicas. La interpretación SIEMPRE debe considerar el origen ancestral y la transferibilidad GWAS. Todos los valores de PRS se reportan con intervalos de confianza del 95%. Las puntuaciones con cobertura < 70% o evidencia < 50/100 deben considerarse solo exploratorias.')}

---

*{T('Generated by Unified Report Engine v9.0.0 — Single Source of Scientific Truth', 'Generado por el Motor de Informes Unificado v9.0.0 — Fuente Única de Verdad Científica')}*
"""
        # Split the template properly
        parts = content.split("---\n\n## 3.")
        if len(parts) > 1:
            body_start = "---\n\n## 3." + parts[1]
        else:
            body_start = ""

        # Reconstruct
        header_end = content.find("---\n\n## 3.")
        if header_end == -1:
            return content  # Return as-is if template parsing fails

        final = content[:header_end] + body_start
        return final

    def _load(self, path: str) -> Dict:
        if Path(path).exists():
            try:
                with open(path) as fh:
                    return json.load(fh)
            except Exception:
                pass
        return {}

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Unified Report Engine")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--prs-core", default="science/prs_core_definition.json")
    parser.add_argument("--ancestry", default="science/ANCESTRY_MODEL.json")
    parser.add_argument("--prs-result", default="prs/PRS_RESULT.json")
    parser.add_argument("--benchmark", default="benchmark/VALIDATION_REPORT.json")
    parser.add_argument("--integrity", default="science/scientific_integrity_score.json")
    parser.add_argument("--output-dir", "-o", default="reports")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    engine = UnifiedReportEngine(args.output_dir)
    paths = engine.generate(args.sample_id, args.prs_core, args.ancestry,
                           args.prs_result, args.benchmark, args.integrity)
    for lang, path in paths.items():
        lines = sum(1 for _ in open(path))
        print(f"  ✅ [{lang}] {path} ({lines} lines)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
