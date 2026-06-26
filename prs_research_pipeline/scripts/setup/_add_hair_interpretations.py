#!/usr/bin/env python3
"""Add bilingual interpretations for new hair color traits."""
import json
from pathlib import Path

BASE = Path(__file__).parent.parent.parent  # prs_research_pipeline

# ── Bilingual hair color interpretations ──
HAIR_INTERPRETATIONS = {
    "Hair color (red)": {
        "en": {
            "trait": "Hair color (red)",
            "name": "Red Hair",
            "gene": "MC1R",
            "pathway": "Melanocortin 1 receptor — melanin type switching (eumelanin ↔ pheomelanin)",
            "description": "Variants in MC1R shift melanin production from dark eumelanin to red/yellow pheomelanin. Having multiple MC1R red hair variants strongly predicts the red hair phenotype. This trait uses 3 MC1R SNPs: rs1805007 (R151C), rs11547464 (R160W), rs1110400 (I155T).",
            "risk_category": "medium",
            "risk_label": "Average",
            "z_score_population": 1.68,
            "percentile_population": 95.3,
            "narrative": "Your genetic profile shows elevated MC1R red hair variant burden relative to European reference distributions. This does NOT mean you have red hair — it indicates you carry genetic variants associated with the red hair phenotype, which may affect sun sensitivity and anesthesia response even without visible red hair.",
            "dietary_context": [
                "MC1R variants reduce tanning ability — use broad-spectrum sunscreen (SPF 50+)",
                "Red hair genetics are linked to altered pain perception and anesthesia requirements",
                "Increased vitamin D synthesis efficiency in low-UV environments — evolutionary advantage"
            ],
            "evidence": "GWAS & candidate gene studies (Valverde 1995, PMID: 7581459; Han 2008, PMID: 18483556)",
            "population_context": "Calibrated against EUR (European) reference. MC1R red hair variants are most common in Northern/Western European populations (~20-30% carry at least one variant)."
        },
        "es": {
            "trait": "Hair color (red)",
            "name": "Pelo Rojo",
            "gene": "MC1R",
            "pathway": "Receptor de melanocortina 1 — cambio de tipo de melanina (eumelanina ↔ feomelanina)",
            "description": "Las variantes en MC1R desplazan la producción de melanina de eumelanina oscura a feomelanina roja/amarilla. Tener múltiples variantes MC1R predice fuertemente el fenotipo pelirrojo. Este rasgo usa 3 SNPs de MC1R: rs1805007 (R151C), rs11547464 (R160W), rs1110400 (I155T).",
            "risk_category": "medium",
            "risk_label": "Promedio",
            "z_score_population": 1.68,
            "percentile_population": 95.3,
            "narrative": "Tu perfil genético muestra una carga elevada de variantes MC1R de pelo rojo en relación con las distribuciones de referencia europeas. Esto NO significa que tengas el pelo rojo — indica que portas variantes genéticas asociadas con el fenotipo pelirrojo, lo que puede afectar la sensibilidad al sol y la respuesta a la anestesia incluso sin pelo rojo visible.",
            "dietary_context": [
                "Las variantes MC1R reducen la capacidad de bronceado — usa protector solar de amplio espectro (FPS 50+)",
                "La genética del pelo rojo está vinculada a una percepción alterada del dolor y requisitos de anestesia",
                "Mayor eficiencia de síntesis de vitamina D en entornos con poca luz UV — ventaja evolutiva"
            ],
            "evidence": "GWAS y estudios de genes candidatos (Valverde 1995, PMID: 7581459; Han 2008, PMID: 18483556)",
            "population_context": "Calibrado contra la referencia EUR (europea). Las variantes MC1R de pelo rojo son más comunes en poblaciones del norte/oeste de Europa (~20-30% portan al menos una variante)."
        }
    },
    "Hair color (blonde)": {
        "en": {
            "trait": "Hair color (blonde)",
            "name": "Blonde Hair",
            "gene": "KITLG, TYR, IRF4, HERC2",
            "pathway": "Melanogenesis — melanin synthesis, transport, and regulation",
            "description": "Blonde hair results from reduced eumelanin production regulated by multiple genes: KITLG (melanocyte development), TYR (tyrosinase enzyme), IRF4 (transcription factor), and HERC2/OCA2 (melanosomal pH regulation). European-specific alleles at these loci contribute to the blonde hair phenotype.",
            "risk_category": "medium",
            "risk_label": "Average",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Your polygenic profile for blonde hair is assessed relative to European reference distributions. Hair color is a cosmetic trait, not a health risk — these results are for ancestry and phenotype curiosity only.",
            "dietary_context": [],
            "evidence": "GWAS (Sulem 2007, PMID: 17999355; Han 2008, PMID: 18483556; Adhikari 2019, PMID: 30895295)",
            "population_context": "Calibrated against EUR (European) reference. Blonde hair is most common in Northern European populations (~5-20% frequency depending on region)."
        },
        "es": {
            "trait": "Hair color (blonde)",
            "name": "Pelo Rubio",
            "gene": "KITLG, TYR, IRF4, HERC2",
            "pathway": "Melanogénesis — síntesis, transporte y regulación de melanina",
            "description": "El pelo rubio resulta de una producción reducida de eumelanina regulada por múltiples genes: KITLG (desarrollo de melanocitos), TYR (enzima tirosinasa), IRF4 (factor de transcripción) y HERC2/OCA2 (regulación del pH melanosomal). Alelos específicos europeos contribuyen al fenotipo rubio.",
            "risk_category": "medium",
            "risk_label": "Promedio",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Tu perfil poligénico para pelo rubio se evalúa en relación con las distribuciones de referencia europeas. El color de pelo es un rasgo cosmético, no un riesgo para la salud — estos resultados son solo por curiosidad sobre ascendencia y fenotipo.",
            "dietary_context": [],
            "evidence": "GWAS (Sulem 2007, PMID: 17999355; Han 2008, PMID: 18483556; Adhikari 2019, PMID: 30895295)",
            "population_context": "Calibrado contra la referencia EUR (europea). El pelo rubio es más común en poblaciones del norte de Europa (~5-20% de frecuencia según la región)."
        }
    },
    "Hair color (brown)": {
        "en": {
            "trait": "Hair color (brown)",
            "name": "Brown Hair",
            "gene": "HERC2, SLC45A2, IRF4",
            "pathway": "Melanogenesis — eumelanin production and distribution",
            "description": "Brown hair is the most common hair color globally and in European populations. It is regulated by HERC2/OCA2 (melanosomal pH), SLC45A2 (melanosomal transporter), and IRF4 (melanocyte transcription factor). These variants influence the amount of dark eumelanin deposited in hair shafts.",
            "risk_category": "medium",
            "risk_label": "Average",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Your polygenic profile for brown hair pigmentation is assessed relative to European reference distributions. This is a cosmetic trait — no health implications.",
            "dietary_context": [],
            "evidence": "GWAS (Han 2008, PMID: 18483556; Sulem 2007, PMID: 17999355)",
            "population_context": "Calibrated against EUR (European) reference. Brown hair is the most common hair color in Southern and Central Europe (60-80% frequency)."
        },
        "es": {
            "trait": "Hair color (brown)",
            "name": "Pelo Castaño",
            "gene": "HERC2, SLC45A2, IRF4",
            "pathway": "Melanogénesis — producción y distribución de eumelanina",
            "description": "El pelo castaño es el color de pelo más común globalmente y en poblaciones europeas. Está regulado por HERC2/OCA2 (pH melanosomal), SLC45A2 (transportador melanosomal) e IRF4 (factor de transcripción melanocítico). Estas variantes influyen en la cantidad de eumelanina oscura depositada en el cabello.",
            "risk_category": "medium",
            "risk_label": "Promedio",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Tu perfil poligénico para pigmentación de pelo castaño se evalúa en relación con distribuciones de referencia europeas. Es un rasgo cosmético — sin implicaciones para la salud.",
            "dietary_context": [],
            "evidence": "GWAS (Han 2008, PMID: 18483556; Sulem 2007, PMID: 17999355)",
            "population_context": "Calibrado contra la referencia EUR (europea). El pelo castaño es el color más común en el sur y centro de Europa (60-80% de frecuencia)."
        }
    },
    "Hair color (black)": {
        "en": {
            "trait": "Hair color (black)",
            "name": "Black Hair",
            "gene": "SLC45A2, SLC24A5",
            "pathway": "Melanogenesis — maximal eumelanin production",
            "description": "Black hair is the most common hair color worldwide, predominant in African, East Asian, South Asian, and Native American populations. It is associated with ancestral (non-derived) alleles at SLC45A2 and SLC24A5, which maintain full eumelanin production. The derived light-skin alleles at these genes swept to near-fixation in European populations but are rare globally.",
            "risk_category": "medium",
            "risk_label": "Average",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Your polygenic profile for black hair pigmentation is assessed using the key pigmentation genes SLC45A2 and SLC24A5. Black hair is the ancestral human trait — lighter hair colors evolved later in specific populations. This is a cosmetic trait with no health implications.",
            "dietary_context": [],
            "evidence": "GWAS (Lamason 2005, PMID: 16357253; Han 2008, PMID: 18483556; Basu Mallick 2013, PMID: 24045820)",
            "population_context": "Calibrated against global reference populations. Black hair is the most common hair color worldwide (>70% of global population)."
        },
        "es": {
            "trait": "Hair color (black)",
            "name": "Pelo Negro",
            "gene": "SLC45A2, SLC24A5",
            "pathway": "Melanogénesis — producción máxima de eumelanina",
            "description": "El pelo negro es el color de pelo más común en el mundo, predominante en poblaciones africanas, asiáticas orientales, del sur de Asia y nativas americanas. Está asociado con alelos ancestrales en SLC45A2 y SLC24A5, que mantienen la producción completa de eumelanina. Los alelos derivados de piel clara en estos genes llegaron a casi fijarse en poblaciones europeas pero son raros globalmente.",
            "risk_category": "medium",
            "risk_label": "Promedio",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Tu perfil poligénico para pigmentación de pelo negro se evalúa usando los genes clave de pigmentación SLC45A2 y SLC24A5. El pelo negro es el rasgo humano ancestral — los colores de pelo más claros evolucionaron después en poblaciones específicas. Es un rasgo cosmético sin implicaciones para la salud.",
            "dietary_context": [],
            "evidence": "GWAS (Lamason 2005, PMID: 16357253; Han 2008, PMID: 18483556; Basu Mallick 2013, PMID: 24045820)",
            "population_context": "Calibrado contra poblaciones de referencia globales. El pelo negro es el color de pelo más común en el mundo (>70% de la población global)."
        }
    },
    "Hair graying (early)": {
        "en": {
            "trait": "Hair graying (early)",
            "name": "Early Hair Graying",
            "gene": "IRF4",
            "pathway": "Melanocyte stem cell maintenance — IRF4 regulates MITF via TFAP2A",
            "description": "IRF4 is the strongest known genetic predictor of early hair graying. The T allele at rs12203592 reduces IRF4 expression in melanocytes, accelerating the depletion of melanocyte stem cells in hair follicles. Each copy of the T allele approximately doubles the odds of early graying (before age 30 in Europeans).",
            "risk_category": "medium",
            "risk_label": "Average",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Your genetic profile for IRF4 is assessed relative to European reference distributions. Early graying is influenced by genetics (~30% heritability) but also strongly by oxidative stress, smoking, vitamin B12 deficiency, and autoimmune conditions. Genetics sets predisposition — lifestyle modulates onset.",
            "dietary_context": [
                "Adequate B12, folate, and copper intake supports melanocyte function",
                "Smoking accelerates hair graying — oxidative stress depletes melanocyte stem cells",
                "Chronic stress has been linked to premature graying in multiple studies"
            ],
            "evidence": "GWAS (Adhikari 2016, PMID: 27182965; largest GWAS of hair graying to date, n=6,357 Latin Americans)",
            "population_context": "Calibrated against EUR (European) reference. The rs12203592 T allele frequency is ~17% in Europeans, <1% in Africans and East Asians."
        },
        "es": {
            "trait": "Hair graying (early)",
            "name": "Canas Prematuras",
            "gene": "IRF4",
            "pathway": "Mantenimiento de células madre melanocíticas — IRF4 regula MITF vía TFAP2A",
            "description": "IRF4 es el predictor genético más fuerte conocido de canas prematuras. El alelo T en rs12203592 reduce la expresión de IRF4 en melanocitos, acelerando el agotamiento de células madre melanocíticas en los folículos pilosos. Cada copia del alelo T aproximadamente duplica la probabilidad de encanecer temprano (antes de los 30 años en europeos).",
            "risk_category": "medium",
            "risk_label": "Promedio",
            "z_score_population": 0.0,
            "percentile_population": 50.0,
            "narrative": "Tu perfil genético para IRF4 se evalúa en relación con distribuciones de referencia europeas. Las canas prematuras están influenciadas por la genética (~30% heredabilidad) pero también fuertemente por estrés oxidativo, tabaquismo, deficiencia de B12 y condiciones autoinmunes. La genética establece la predisposición — el estilo de vida modula el inicio.",
            "dietary_context": [
                "Una ingesta adecuada de B12, folato y cobre apoya la función melanocítica",
                "Fumar acelera las canas — el estrés oxidativo agota las células madre melanocíticas",
                "El estrés crónico se ha vinculado a canas prematuras en múltiples estudios"
            ],
            "evidence": "GWAS (Adhikari 2016, PMID: 27182965; el GWAS más grande de canas hasta la fecha, n=6,357 latinoamericanos)",
            "population_context": "Calibrado contra la referencia EUR (europea). La frecuencia del alelo T de rs12203592 es ~17% en europeos, <1% en africanos y asiáticos orientales."
        }
    }
}

# ── Update both files ──
for lang_file, default_lang in [("interpretations_en.json", "en"), ("interpretations_es.json", "es")]:
    path = BASE / "interpretations" / lang_file
    data = json.loads(path.read_text())

    for trait, texts in HAIR_INTERPRETATIONS.items():
        if default_lang == "en":
            # For EN file, use the EN text directly
            data["trait_interpretations"][trait] = texts["en"]
        else:
            # For ES file, use the ES text
            if trait in texts:
                data["trait_interpretations"][trait] = texts["es"]
            else:
                # Fallback to EN for any missing ES
                data["trait_interpretations"][trait] = texts["en"]

    path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print(f"Updated {lang_file}: {len(data['trait_interpretations'])} traits")

print("\nNew hair traits added to interpretations:")
for trait in HAIR_INTERPRETATIONS:
    print(f"  • {trait}")
