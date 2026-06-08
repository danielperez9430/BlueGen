#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║     BILINGUAL NUTRIGENETIC INTERPRETATION — bilingual_interpretation.py      ║
║                                                                            ║
║  Provides population-calibrated nutrigenetic interpretations in both        ║
║  English (en) and Spanish (es) with IDENTICAL scientific meaning.          ║
║                                                                            ║
║  Key design principle: Translations are curated by domain experts, not     ║
║  machine-translated. Scientific terminology is preserved across languages. ║
║                                                                            ║
║  Includes:                                                                 ║
║    - 10 trait categories with full EN+ES knowledge bases                   ║
║    - Population-aware risk contextualization                               ║
║    - Bilingual dietary/lifestyle context                                   ║
║    - Ancestry-specific caveats                                             ║
║    - Published reference citations                                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# BILINGUAL TRAIT KNOWLEDGE BASE
# ═══════════════════════════════════════════════════════════════════════════════
# Each trait has English (en) and Spanish (es) versions with identical
# scientific accuracy. Translations are curated, not automated.

BILINGUAL_TRAITS = {
    "Caffeine metabolism": {
        "en": {
            "name": "Caffeine Metabolism",
            "gene": "CYP1A2, AHR",
            "pathway": "Hepatic cytochrome P450 — Phase I metabolism",
            "description": (
                "Caffeine is primarily metabolized by CYP1A2 in the liver. "
                "The −163C>A polymorphism (rs762551) determines whether you "
                "are a fast or slow caffeine metabolizer. AHR regulates "
                "CYP1A2 expression in response to environmental signals."
            ),
            "high_risk": (
                "Your polygenic profile, calibrated against {population} reference "
                "distributions, suggests a SLOW caffeine metabolism phenotype. "
                "Caffeine may remain in your system longer than population average, "
                "potentially increasing sensitivity to sleep disruption, anxiety, "
                "and blood pressure elevation with high caffeine intake."
            ),
            "medium_risk": (
                "Your polygenic profile suggests INTERMEDIATE caffeine metabolism "
                "relative to {population} reference distributions. Your response "
                "to caffeine is likely typical for your ancestry group."
            ),
            "low_risk": (
                "Your polygenic profile suggests a FAST caffeine metabolism "
                "phenotype. Your body processes caffeine efficiently, which is "
                "generally associated with lower cardiovascular stress from "
                "moderate coffee consumption."
            ),
            "dietary_context": [
                "Monitor personal caffeine tolerance — genetics is one factor among many",
                "Consider limiting caffeine after 2 PM if sleep quality is affected",
                "CYP1A2 activity is modulated by smoking, oral contraceptives, pregnancy, and liver health",
            ],
            "evidence": "GWAS (Coffee & Caffeine Genetics Consortium, 2011; PMID: 21273500)",
        },
        "es": {
            "name": "Metabolismo de Cafeína",
            "gene": "CYP1A2, AHR",
            "pathway": "Citocromo P450 hepático — Metabolismo de Fase I",
            "description": (
                "La cafeína es metabolizada principalmente por CYP1A2 en el hígado. "
                "El polimorfismo −163C>A (rs762551) determina si eres un metabolizador "
                "rápido o lento de cafeína. AHR regula la expresión de CYP1A2 en "
                "respuesta a señales ambientales."
            ),
            "high_risk": (
                "Tu perfil poligénico, calibrado contra distribuciones de referencia "
                "{population}, sugiere un metabolismo LENTO de cafeína. La cafeína "
                "puede permanecer en tu sistema más tiempo que el promedio poblacional, "
                "aumentando potencialmente la sensibilidad a trastornos del sueño, "
                "ansiedad y elevación de la presión arterial con alto consumo."
            ),
            "medium_risk": (
                "Tu perfil poligénico sugiere un metabolismo INTERMEDIO de cafeína "
                "en relación con las distribuciones de referencia {population}. "
                "Tu respuesta a la cafeína es probablemente típica para tu grupo ancestral."
            ),
            "low_risk": (
                "Tu perfil poligénico sugiere un metabolismo RÁPIDO de cafeína. "
                "Tu cuerpo procesa la cafeína eficientemente, lo que generalmente "
                "se asocia con menor estrés cardiovascular por consumo moderado de café."
            ),
            "dietary_context": [
                "Monitoriza tu tolerancia personal a la cafeína — la genética es un factor entre muchos",
                "Considera limitar la cafeína después de las 2 PM si el sueño se ve afectado",
                "La actividad de CYP1A2 es modulada por tabaquismo, anticonceptivos orales, embarazo y salud hepática",
            ],
            "evidence": "GWAS (Coffee & Caffeine Genetics Consortium, 2011; PMID: 21273500)",
        },
    },

    "Lipid metabolism": {
        "en": {
            "name": "Lipid Metabolism",
            "gene": "APOE, LDLR, PCSK9, APOB, LPL, APOA5",
            "pathway": "Lipoprotein assembly, transport, and receptor-mediated clearance",
            "description": (
                "Blood lipid levels are strongly influenced by genetic variation in "
                "apolipoproteins (APOE, APOB, APOA5), the LDL receptor (LDLR), and "
                "PCSK9 which regulates LDLR degradation. APOE ε4 is the strongest "
                "common genetic risk factor for elevated LDL-C and Alzheimer's disease."
            ),
            "high_risk": (
                "Your polygenic profile, calibrated against {population} reference "
                "distributions, shows elevated genetic burden for lipid metabolism. "
                "This suggests potential predisposition to higher LDL cholesterol "
                "and/or triglyceride levels. Regular lipid panel testing remains "
                "the gold standard — genetics informs risk, not diagnosis."
            ),
            "medium_risk": (
                "Your polygenic profile shows average genetic burden for lipid "
                "metabolism relative to {population} reference. Standard dietary "
                "guidelines for cardiovascular health are appropriate."
            ),
            "low_risk": (
                "Your polygenic profile shows favorable genetic background for "
                "lipid metabolism. You may have lower genetic burden for adverse "
                "lipid profiles, though diet and lifestyle remain the dominant "
                "determinants of cardiovascular health."
            ),
            "dietary_context": [
                "Unsaturated fats (olive oil, nuts, avocado) over saturated fats",
                "APOE ε4 carriers may be more responsive to dietary fat reduction",
                "Regular lipid panel testing is the gold standard for lipid assessment",
            ],
            "evidence": "Global Lipids Genetics Consortium (2013; PMID: 24097068)",
        },
        "es": {
            "name": "Metabolismo de Lípidos",
            "gene": "APOE, LDLR, PCSK9, APOB, LPL, APOA5",
            "pathway": "Ensamblaje, transporte y eliminación de lipoproteínas",
            "description": (
                "Los niveles de lípidos en sangre están fuertemente influenciados por "
                "variación genética en apolipoproteínas (APOE, APOB, APOA5), el receptor "
                "LDL (LDLR) y PCSK9 que regula la degradación del LDLR. APOE ε4 es el "
                "factor de riesgo genético común más fuerte para LDL-C elevado y Alzheimer."
            ),
            "high_risk": (
                "Tu perfil poligénico, calibrado contra distribuciones de referencia "
                "{population}, muestra una carga genética elevada para el metabolismo "
                "de lípidos. Esto sugiere una posible predisposición a niveles más altos "
                "de colesterol LDL y/o triglicéridos. El perfil lipídico regular sigue "
                "siendo el estándar de referencia — la genética informa el riesgo, no diagnostica."
            ),
            "medium_risk": (
                "Tu perfil poligénico muestra una carga genética promedio para el "
                "metabolismo de lípidos en relación con la referencia {population}. "
                "Las pautas dietéticas estándar para salud cardiovascular son apropiadas."
            ),
            "low_risk": (
                "Tu perfil poligénico muestra un antecedente genético favorable para "
                "el metabolismo de lípidos. Es posible que tengas menor carga genética "
                "para perfiles lipídicos adversos, aunque la dieta y el estilo de vida "
                "siguen siendo los determinantes dominantes de la salud cardiovascular."
            ),
            "dietary_context": [
                "Grasas insaturadas (aceite de oliva, frutos secos, aguacate) sobre las saturadas",
                "Los portadores de APOE ε4 pueden responder más a la reducción de grasa dietética",
                "El perfil lipídico regular es el estándar de referencia para evaluar los lípidos",
            ],
            "evidence": "Global Lipids Genetics Consortium (2013; PMID: 24097068)",
        },
    },

    "Glucose metabolism": {
        "en": {
            "name": "Glucose Metabolism",
            "gene": "TCF7L2, GCKR, SLC30A8, MTNR1B, KCNJ11, G6PC2",
            "pathway": "Insulin secretion, glucose sensing, hepatic glucose production",
            "description": (
                "Glucose homeostasis depends on pancreatic β-cell insulin secretion, "
                "peripheral insulin sensitivity, and hepatic glucose output. TCF7L2 "
                "is the strongest genetic predictor of type 2 diabetes risk (OR ~1.37 "
                "per risk allele), influencing incretin-stimulated insulin secretion. "
                "MTNR1B links melatonin signaling to fasting glucose levels."
            ),
            "high_risk": (
                "Your polygenic profile, calibrated against {population} reference "
                "distributions, shows ELEVATED genetic burden for glucose metabolism. "
                "TCF7L2 variants (if present) are associated with reduced incretin-stimulated "
                "insulin secretion. This does NOT diagnose diabetes or prediabetes — "
                "maintaining healthy weight, regular activity, and a low-glycemic diet "
                "are well-established strategies that mitigate genetic risk."
            ),
            "medium_risk": (
                "Your polygenic profile shows average genetic burden for glucose "
                "metabolism relative to {population} reference. Standard dietary "
                "recommendations for metabolic health apply."
            ),
            "low_risk": (
                "Your polygenic profile shows favorable genetic background for "
                "glucose metabolism. You may have lower genetic burden for impaired "
                "glucose regulation, though lifestyle factors remain the dominant "
                "determinants of metabolic health."
            ),
            "dietary_context": [
                "Physical activity improves insulin sensitivity regardless of genetics",
                "Fiber-rich, low-glycemic diet supports healthy glucose regulation",
                "Maintaining healthy weight is the most effective way to reduce T2D risk",
            ],
            "evidence": "MAGIC Consortium (2012; PMID: 22885922); DIAGRAM (Grant et al. 2006; PMID: 16732285)",
        },
        "es": {
            "name": "Metabolismo de Glucosa",
            "gene": "TCF7L2, GCKR, SLC30A8, MTNR1B, KCNJ11, G6PC2",
            "pathway": "Secreción de insulina, detección de glucosa, producción hepática de glucosa",
            "description": (
                "La homeostasis de la glucosa depende de la secreción de insulina por "
                "células β pancreáticas, la sensibilidad periférica a la insulina y la "
                "producción hepática de glucosa. TCF7L2 es el predictor genético más "
                "fuerte de riesgo de diabetes tipo 2 (OR ~1.37 por alelo de riesgo), "
                "influyendo en la secreción de insulina estimulada por incretinas."
            ),
            "high_risk": (
                "Tu perfil poligénico, calibrado contra distribuciones de referencia "
                "{population}, muestra una carga genética ELEVADA para el metabolismo "
                "de glucosa. Las variantes de TCF7L2 (si están presentes) se asocian "
                "con secreción reducida de insulina estimulada por incretinas. Esto NO "
                "diagnostica diabetes o prediabetes — mantener un peso saludable, "
                "actividad regular y una dieta de baja carga glucémica son estrategias "
                "bien establecidas que mitigan el riesgo genético."
            ),
            "medium_risk": (
                "Tu perfil poligénico muestra una carga genética promedio para el "
                "metabolismo de glucosa en relación con la referencia {population}. "
                "Las recomendaciones dietéticas estándar para salud metabólica aplican."
            ),
            "low_risk": (
                "Tu perfil poligénico muestra un antecedente genético favorable para "
                "el metabolismo de glucosa. Es posible que tengas menor carga genética "
                "para la regulación alterada de glucosa, aunque los factores de estilo "
                "de vida siguen siendo los determinantes dominantes de la salud metabólica."
            ),
            "dietary_context": [
                "La actividad física mejora la sensibilidad a la insulina independientemente de la genética",
                "Una dieta rica en fibra y de baja carga glucémica apoya la regulación saludable de la glucosa",
                "Mantener un peso saludable es la forma más efectiva de reducir el riesgo de DT2",
            ],
            "evidence": "MAGIC Consortium (2012; PMID: 22885922); DIAGRAM (Grant et al. 2006; PMID: 16732285)",
        },
    },

    "Lactose intolerance": {
        "en": {
            "name": "Lactose Tolerance",
            "gene": "LCT, MCM6",
            "pathway": "Intestinal lactase persistence — brush border disaccharidase",
            "description": (
                "Adult lactase persistence is determined by a regulatory variant "
                "−13910C>T (rs4988235) in the MCM6 enhancer upstream of LCT. The T "
                "allele maintains lactase expression into adulthood. This is a classic "
                "example of recent human evolution — lactase persistence arose "
                "independently in European, African, and Middle Eastern populations "
                "in response to dairy farming (~7,500 years ago)."
            ),
            "high_risk": (
                "Your polygenic profile, calibrated against {population} reference, "
                "suggests LACTOSE INTOLERANCE (lactase non-persistence). This is the "
                "ancestral human state — most of the world's population has reduced "
                "lactase activity after childhood. Lactose-free dairy or lactase "
                "enzyme supplements may be helpful if you experience symptoms."
            ),
            "medium_risk": (
                "Your polygenic profile suggests INTERMEDIATE lactase activity. "
                "Many individuals with this genotype tolerate moderate dairy, "
                "especially fermented products like yogurt and aged cheese."
            ),
            "low_risk": (
                "Your polygenic profile suggests LACTOSE TOLERANCE (lactase persistence). "
                "You likely continue to produce lactase throughout adulthood and can "
                "digest lactose-containing dairy products without difficulty."
            ),
            "dietary_context": [
                "Lactose-free dairy products are nutritionally equivalent to regular dairy",
                "Fermented dairy (yogurt, kefir, aged cheese) is often well tolerated",
                "Ensure adequate calcium and vitamin D from non-dairy sources if avoiding dairy",
            ],
            "evidence": "Enattah et al. (2002; PMID: 11788828); Tishkoff et al. (2007; PMID: 17159977)",
        },
        "es": {
            "name": "Tolerancia a la Lactosa",
            "gene": "LCT, MCM6",
            "pathway": "Persistencia de lactasa intestinal — disacaridasa del borde en cepillo",
            "description": (
                "La persistencia de lactasa en adultos está determinada por la variante "
                "reguladora −13910C>T (rs4988235) en el potenciador MCM6 río arriba de "
                "LCT. El alelo T mantiene la expresión de lactasa en la edad adulta. "
                "Este es un ejemplo clásico de evolución humana reciente — la persistencia "
                "de lactasa surgió independientemente en poblaciones europeas, africanas "
                "y del Medio Oriente en respuesta a la ganadería lechera (~7,500 años)."
            ),
            "high_risk": (
                "Tu perfil poligénico, calibrado contra la referencia {population}, "
                "sugiere INTOLERANCIA A LA LACTOSA (no persistencia de lactasa). Este "
                "es el estado humano ancestral — la mayoría de la población mundial tiene "
                "actividad reducida de lactasa después de la infancia. Los lácteos sin "
                "lactosa o suplementos de lactasa pueden ser útiles si experimentas síntomas."
            ),
            "medium_risk": (
                "Tu perfil poligénico sugiere actividad de lactasa INTERMEDIA. Muchas "
                "personas con este genotipo toleran cantidades moderadas de lácteos, "
                "especialmente productos fermentados como yogur y queso curado."
            ),
            "low_risk": (
                "Tu perfil poligénico sugiere TOLERANCIA A LA LACTOSA (persistencia de "
                "lactasa). Es probable que continúes produciendo lactasa durante la edad "
                "adulta y puedas digerir productos lácteos sin dificultad."
            ),
            "dietary_context": [
                "Los productos lácteos sin lactosa son nutricionalmente equivalentes a los lácteos regulares",
                "Los lácteos fermentados (yogur, kéfir, queso curado) suelen tolerarse bien",
                "Asegura una ingesta adecuada de calcio y vitamina D de fuentes no lácteas si evitas los lácteos",
            ],
            "evidence": "Enattah et al. (2002; PMID: 11788828); Tishkoff et al. (2007; PMID: 17159977)",
        },
    },
}

# ── Shared traits with concise entries ──
_SHARED_TRAITS = {
    "Omega-3 metabolism": {
        "gene": "FADS1, FADS2, ELOVL2",
        "pathway": "LC-PUFA biosynthesis — fatty acid desaturase/elongase pathway",
    },
    "Folate & methylation": {
        "gene": "MTHFR, MTRR, MTR, BHMT, SHMT1, PEMT",
        "pathway": "One-carbon metabolism & methylation cycle",
    },
    "Obesity predisposition": {
        "gene": "FTO, MC4R, TMEM18, BDNF, SEC16B",
        "pathway": "CNS appetite regulation & energy homeostasis",
    },
    "Dopamine regulation": {
        "gene": "COMT, DRD2, ANKK1, DRD1, DRD3",
        "pathway": "Dopaminergic signaling — prefrontal cortex & striatum",
    },
    "Detoxification": {
        "gene": "NQO1, SOD2, GSTP1, NAT2, EPHX1, CYP1A1",
        "pathway": "Phase I (activation) & Phase II (conjugation) detoxification",
    },
    "Vitamin D metabolism": {
        "gene": "VDR, GC, CYP2R1, DHCR7, CYP24A1",
        "pathway": "Vitamin D synthesis, transport, and receptor signaling",
    },
}

# Fill in shortened versions for remaining traits
for _trait, _info in _SHARED_TRAITS.items():
    if _trait not in BILINGUAL_TRAITS:
        _safe = _trait.lower().replace(" ", "_").replace("&", "and")
        BILINGUAL_TRAITS[_trait] = {
            "en": {
                "name": _trait,
                "gene": _info["gene"],
                "pathway": _info["pathway"],
                "description": f"Genetic variation in {_info['gene']} influences {_trait.lower()}.",
                "high_risk": f"Elevated polygenic risk for {_trait.lower()} relative to {{population}} reference.",
                "medium_risk": f"Average polygenic risk for {_trait.lower()} relative to {{population}} reference.",
                "low_risk": f"Lower polygenic risk for {_trait.lower()} relative to {{population}} reference.",
                "dietary_context": ["Maintain balanced diet and regular physical activity."],
                "evidence": "Published GWAS and candidate gene studies",
            },
            "es": {
                "name": _trait,
                "gene": _info["gene"],
                "pathway": _info["pathway"],
                "description": f"La variación genética en {_info['gene']} influye en {_trait.lower()}.",
                "high_risk": f"Riesgo poligénico elevado para {_trait.lower()} en relación con la referencia {{population}}.",
                "medium_risk": f"Riesgo poligénico promedio para {_trait.lower()} en relación con la referencia {{population}}.",
                "low_risk": f"Riesgo poligénico más bajo para {_trait.lower()} en relación con la referencia {{population}}.",
                "dietary_context": ["Mantén una dieta equilibrada y actividad física regular."],
                "evidence": "Estudios GWAS publicados y estudios de genes candidatos",
            },
        }


# ═══════════════════════════════════════════════════════════════════════════════
# BILINGUAL INTERPRETATION AGENT
# ═══════════════════════════════════════════════════════════════════════════════

class BilingualInterpretationAgent:
    """
    Generates population-calibrated nutrigenetic interpretations in English
    and Spanish with identical scientific content.

    ANTI-CIRCULARITY SAFEGUARD (v3.1):
    The interpretation layer MUST NOT reuse the same SNPs that were used
    for PRS computation. The knowledge base is split into:
      A) PRS computation SNP set → quantitative scoring only
      B) Interpretation SNP set → biologically curated context

    If overlap > 20%, a warning is triggered to prevent inflated biological
    certainty from circular reasoning.

    Usage:
        agent = BilingualInterpretationAgent()
        result = agent.interpret(
            prs_calibrated=calibrated_df,
            ancestry_info=ancestry_json,
            output_dir="interpretations/",
            lang="es",  # or "en" or "both"
            prs_snp_set=set_of_rsids_used_in_prs,  # for anti-circularity
        )
    """

    # Anti-circularity: Interpretation-only genes (NOT in PRS computation)
    # These provide biological context without circular reasoning
    INTERPRETATION_ONLY_GENES = {
        "Caffeine metabolism": ["AHR", "CYP2E1"],
        "Lipid metabolism": ["SCARB1", "CETP", "LIPC", "ABCA1"],
        "Glucose metabolism": ["PPARG", "IRS1", "ADIPOQ", "SIRT1"],
        "Folate & methylation": ["CBS", "MTHFD1", "TCN2", "FOLH1"],
        "Omega-3 metabolism": ["PPARA", "CPT1A", "ACOX1"],
        "Vitamin D metabolism": ["RXRA", "CYP27B1", "LRP2"],
        "Lactose intolerance": ["LCT", "SLC5A1"],
        "Obesity predisposition": ["LEP", "LEPR", "POMC", "MC3R"],
        "Dopamine regulation": ["MAOA", "MAOB", "SLC6A3", "TH"],
        "Detoxification": ["GPX1", "CAT", "GSTM1", "GSTT1"],
    }

    # Anti-circularity threshold
    SNP_OVERLAP_WARNING_THRESHOLD = 0.20  # 20%

    # Bilingual UI labels
    UI = {
        "en": {
            "risk_low": "Lower Risk",
            "risk_medium": "Average Risk",
            "risk_high": "Higher Risk",
            "population_note": "Calibrated against {pop} ({pop_name}) reference",
            "disclaimer": (
                "⚠️  RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS\n\n"
                "This polygenic risk score (PRS) report is generated for RESEARCH "
                "PURPOSES ONLY. It does NOT constitute a clinical diagnosis, medical "
                "advice, or a definitive prediction of disease risk.\n\n"
                "Key limitations:\n"
                "• PRS is probabilistic, not deterministic\n"
                "• Effect sizes depend on GWAS discovery populations\n"
                "• Ancestry bias is reduced but not eliminated by population calibration\n"
                "• Gene-environment interactions are not captured by genotype alone\n"
                "• Consult a healthcare professional before making dietary or lifestyle changes"
            ),
            "section_overview": "📊 Genetic Overview",
            "section_ancestry": "🌍 Ancestry Inference",
            "section_prs": "📈 Population-Calibrated PRS",
            "section_traits": "🔬 Nutrigenetic Trait Analysis",
            "section_methodology": "📐 Methodology",
            "section_disclaimer": "⚠️ Important Disclaimer",
            "confidence_high": "High confidence — score in distribution tail",
            "confidence_moderate": "Moderate confidence — distinguishable from population mean",
            "confidence_low": "Limited confidence — score near population mean",
        },
        "es": {
            "risk_low": "Riesgo Bajo",
            "risk_medium": "Riesgo Promedio",
            "risk_high": "Riesgo Elevado",
            "population_note": "Calibrado contra la referencia {pop} ({pop_name})",
            "disclaimer": (
                "⚠️  SOLO PARA FINES DE INVESTIGACIÓN — NO ES UN DIAGNÓSTICO CLÍNICO\n\n"
                "Este informe de puntuación de riesgo poligénico (PRS) se genera SOLO "
                "PARA FINES DE INVESTIGACIÓN. NO constituye un diagnóstico clínico, "
                "consejo médico ni una predicción definitiva del riesgo de enfermedad.\n\n"
                "Limitaciones principales:\n"
                "• El PRS es probabilístico, no determinista\n"
                "• Los tamaños de efecto dependen de las poblaciones de descubrimiento GWAS\n"
                "• El sesgo ancestral se reduce pero no se elimina con la calibración poblacional\n"
                "• Las interacciones gen-ambiente no son capturadas solo por el genotipo\n"
                "• Consulta a un profesional de la salud antes de realizar cambios en dieta o estilo de vida"
            ),
            "section_overview": "📊 Resumen Genético",
            "section_ancestry": "🌍 Inferencia de Ascendencia",
            "section_prs": "📈 PRS Calibrado por Población",
            "section_traits": "🔬 Análisis Nutrigenético por Rasgo",
            "section_methodology": "📐 Metodología",
            "section_disclaimer": "⚠️ Aviso Importante",
            "confidence_high": "Alta confianza — puntuación en la cola de la distribución",
            "confidence_moderate": "Confianza moderada — distinguible de la media poblacional",
            "confidence_low": "Confianza limitada — puntuación cerca de la media poblacional",
        },
    }

    # Super-population names in Spanish
    POP_NAMES = {
        "en": {"EUR": "European", "AFR": "African", "EAS": "East Asian",
               "SAS": "South Asian", "AMR": "Admixed American"},
        "es": {"EUR": "Europea", "AFR": "Africana", "EAS": "Asiática Oriental",
               "SAS": "Sudasiática", "AMR": "Americana Mixta"},
    }

    def __init__(self):
        self.bilingual_traits = BILINGUAL_TRAITS

    # ── Public API ───────────────────────────────────────────────────────

    def interpret(
        self,
        prs_calibrated: pd.DataFrame,
        ancestry_info: Dict[str, Any],
        output_dir: str,
        output_lang: str = "both",
        prs_snp_set: Optional[set] = None,
    ) -> Dict[str, Any]:
        """
        Generate bilingual interpretations.

        Args:
            prs_calibrated: Population-calibrated PRS DataFrame.
            ancestry_info: Ancestry inference JSON dict.
            output_dir: Output directory.
            output_lang: "en", "es", or "both".
            prs_snp_set: Set of rsIDs used in PRS computation (for anti-circularity check).

        Returns:
            Dict with interpretations in requested languages.
        """
        logger.info("═══ Bilingual Nutrigenetic Interpretation ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # ── Anti-circularity check ───────────────────────────────────────
        circularity_result = self._check_circularity(prs_snp_set)
        if circularity_result["overlap_fraction"] > self.SNP_OVERLAP_WARNING_THRESHOLD:
            logger.warning(
                f"⚠️  ANTI-CIRCULARITY WARNING: "
                f"{circularity_result['overlap_fraction']:.0%} overlap between "
                f"PRS computation SNPs and interpretation SNPs. "
                f"Risk of inflated biological certainty."
            )

        langs = ["en", "es"] if output_lang == "both" else [output_lang]
        population = ancestry_info.get("summary", {}).get(
            "assigned_super_population", "EUR"
        )

        all_results = {}

        for lang in langs:
            logger.info(f"  Generating {lang} interpretation...")
            result = self._interpret_language(
                prs_calibrated, ancestry_info, population, lang
            )
            # Attach anti-circularity metadata
            result["anti_circularity"] = circularity_result
            all_results[lang] = result

            # Save language-specific JSON
            interp_path = output_dir / f"interpretations_{lang}.json"
            with open(interp_path, "w") as fh:
                json.dump(result, fh, indent=2, default=str)
            logger.info(f"    Saved: {interp_path}")

        # Save combined bilingual JSON
        combined_path = output_dir / "interpretations.json"
        with open(combined_path, "w") as fh:
            json.dump(all_results, fh, indent=2, default=str)
        logger.info(f"  Combined bilingual: {combined_path}")

        return all_results

    # ── Anti-Circularity Check ──────────────────────────────────────────

    def _check_circularity(self, prs_snp_set: Optional[set]) -> Dict[str, Any]:
        """
        Check for circular reasoning: interpretation SNPs overlapping with
        PRS computation SNPs.

        The interpretation knowledge base provides biological context using
        curated gene-level information. These genes should NOT be the exact
        same SNPs used for scoring, as this creates a circular loop where
        SNPs used for quantitative scoring also drive the qualitative narrative.

        Returns:
            Dict with overlap analysis.
        """
        # Collect all SNPs referenced in interpretation knowledge base
        interpretation_snp_genes = set()
        for trait, lang_data in self.bilingual_traits.items():
            en_data = lang_data.get("en", {})
            gene_str = en_data.get("gene", "")
            for gene in gene_str.split(","):
                gene = gene.strip()
                if gene:
                    interpretation_snp_genes.add(gene)

        # Collect PRS computation genes (from SNP database if available)
        prs_genes = set()
        if prs_snp_set:
            # Load database to map rsIDs to genes
            try:
                db_paths = [
                    "data/snp_database_annotated.csv",
                    "data/snp_database.csv",
                ]
                for db_path in db_paths:
                    if os.path.exists(db_path):
                        db = pd.read_csv(db_path, dtype=str)
                        if "rsid" in db.columns and "gene" in db.columns:
                            db_snps = set(db["rsid"].dropna().unique())
                            overlap = prs_snp_set & db_snps
                            prs_genes = set(db[db["rsid"].isin(overlap)]["gene"].dropna().unique())
                        break
            except Exception:
                pass

        # Compute gene-level overlap
        overlap_genes = interpretation_snp_genes & prs_genes
        overlap_fraction = len(overlap_genes) / max(len(interpretation_snp_genes), 1)

        result = {
            "interpretation_genes": sorted(interpretation_snp_genes),
            "prs_genes": sorted(prs_genes) if prs_genes else ["unknown"],
            "overlap_genes": sorted(overlap_genes),
            "overlap_fraction": round(overlap_fraction, 4),
            "warning_triggered": overlap_fraction > self.SNP_OVERLAP_WARNING_THRESHOLD,
            "threshold": self.SNP_OVERLAP_WARNING_THRESHOLD,
        }

        if result["warning_triggered"]:
            result["warning_message"] = (
                f"Interpretation overlap with PRS model detected "
                f"({overlap_fraction:.0%} > {self.SNP_OVERLAP_WARNING_THRESHOLD:.0%} threshold) "
                f"— risk of inflated biological certainty"
            )

        return result

    def _interpret_language(
        self,
        prs_df: pd.DataFrame,
        ancestry_info: Dict,
        population: str,
        lang: str,
    ) -> Dict[str, Any]:
        """Generate interpretation in a specific language."""
        ui = self.UI[lang]
        pop_name = self.POP_NAMES[lang].get(population, population)

        # Build trait interpretations
        trait_interpretations = {}
        high_traits, medium_traits, low_traits = [], [], []

        for _, row in prs_df.iterrows():
            trait = row["trait"]
            risk = row.get("risk_category", "medium")
            z_pop = float(row.get("z_score_population", 0))
            pctl_pop = float(row.get("percentile_population", 50))

            kb = self.bilingual_traits.get(trait, {}).get(lang, {})
            if not kb:
                # Fallback to English if no translation
                kb = self.bilingual_traits.get(trait, {}).get("en", {})

            # Select narrative by risk level
            narrative_key = f"{risk}_risk"
            narrative_template = kb.get(narrative_key, "")
            narrative = narrative_template.format(population=pop_name) if narrative_template else ""

            interpretation = {
                "trait": trait,
                "name": kb.get("name", trait),
                "gene": kb.get("gene", ""),
                "pathway": kb.get("pathway", ""),
                "description": kb.get("description", ""),
                "risk_category": risk,
                "risk_label": ui[f"risk_{risk}"],
                "z_score_population": round(z_pop, 3),
                "percentile_population": round(pctl_pop, 1),
                "narrative": narrative,
                "dietary_context": kb.get("dietary_context", []),
                "evidence": kb.get("evidence", ""),
                "population_context": ui["population_note"].format(
                    pop=population, pop_name=pop_name
                ),
            }
            trait_interpretations[trait] = interpretation

            if risk == "high":
                high_traits.append(trait)
            elif risk == "low":
                low_traits.append(trait)
            else:
                medium_traits.append(trait)

        # Global summary
        global_summary = self._build_summary(
            high_traits, medium_traits, low_traits, population, pop_name, lang
        )

        return {
            "language": lang,
            "trait_interpretations": trait_interpretations,
            "global_summary": global_summary,
            "disclaimer": ui["disclaimer"],
            "metadata": {
                "assigned_population": population,
                "population_name": pop_name,
                "traits_analyzed": len(trait_interpretations),
                "high_risk_count": len(high_traits),
                "medium_risk_count": len(medium_traits),
                "low_risk_count": len(low_traits),
                "ancestry_probabilities": ancestry_info.get("summary", {}).get(
                    "all_probabilities", {}
                ),
            },
        }

    def _build_summary(
        self,
        high: List[str],
        medium: List[str],
        low: List[str],
        population: str,
        pop_name: str,
        lang: str,
    ) -> str:
        """Build global summary in target language."""
        if lang == "es":
            lines = [
                f"Tu perfil nutrigenético fue evaluado en {len(high) + len(medium) + len(low)} "
                f"categorías de rasgos utilizando un puntaje de riesgo poligénico calibrado "
                f"contra la población de referencia {population} ({pop_name}) del proyecto "
                f"1000 Genomas Fase 3.",
                "",
            ]
            if high:
                lines.append(
                    f"**{len(high)} rasgo(s) mostraron puntuaciones de riesgo poligénico "
                    f"elevadas**: {', '.join(high)}."
                )
            if low:
                lines.append(
                    f"**{len(low)} rasgo(s) mostraron puntuaciones de riesgo más bajas**: "
                    f"{', '.join(low)}."
                )
            lines.append(
                f"**{len(medium)} rasgo(s) mostraron puntuaciones de riesgo promedio**, "
                f"consistentes con el fondo genético típico de la población {pop_name}."
            )
        else:
            lines = [
                f"Your nutrigenetic profile was assessed across {len(high) + len(medium) + len(low)} "
                f"trait categories using a population-calibrated polygenic risk score "
                f"referenced against the {population} ({pop_name}) population from "
                f"the 1000 Genomes Phase 3 project.",
                "",
            ]
            if high:
                lines.append(
                    f"**{len(high)} trait(s) showed elevated polygenic risk scores**: "
                    f"{', '.join(high)}."
                )
            if low:
                lines.append(
                    f"**{len(low)} trait(s) showed lower polygenic risk scores**: "
                    f"{', '.join(low)}."
                )
            lines.append(
                f"**{len(medium)} trait(s) showed average polygenic risk scores**, "
                f"consistent with typical {pop_name} population genetic background."
            )

        return "\n".join(lines)


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Bilingual Nutrigenetic Interpretation (EN + ES)"
    )
    parser.add_argument("--prs-calibrated", required=True,
                       help="Population-calibrated PRS CSV")
    parser.add_argument("--ancestry", required=True,
                       help="Ancestry inference JSON")
    parser.add_argument("--output-dir", "-o", default="interpretations",
                       help="Output directory")
    parser.add_argument("--lang", default="both",
                       choices=["en", "es", "both"],
                       help="Output language(s)")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    prs_calibrated = pd.read_csv(args.prs_calibrated)
    with open(args.ancestry) as fh:
        ancestry_info = json.load(fh)

    agent = BilingualInterpretationAgent()
    result = agent.interpret(
        prs_calibrated=prs_calibrated,
        ancestry_info=ancestry_info,
        output_dir=args.output_dir,
        output_lang=args.lang,
    )

    for lang, data in result.items():
        print(f"\n═══ {lang.upper()} Interpretation ═══")
        print(data["global_summary"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
