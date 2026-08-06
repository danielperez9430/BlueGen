"""
Gene/disease -> body-system classifier for grouping ClinVar findings
(IMPROVEMENT_PLAN.md 1.3 - "agrupar por sistema/organo").

There is no single authoritative disease->organ-system taxonomy for
ClinVar's full gene list (genes are pleiotropic; MedGen's own semantic
types are too coarse - "Disease or Syndrome" covers everything). This
module is a pragmatic, extensible heuristic, not a clinical-grade
ontology: gene symbol first (most reliable signal for a single ClinVar
row, since one row is usually one gene), falling back to keyword
matching on the MedGen disease name, falling back to "other" if neither
matches. Coverage is intentionally partial - an unclassified row is
labeled "Other / Multisystem" rather than guessed at.
"""

# Ordered so the report always lists systems in the same, clinically
# conventional order regardless of dict iteration order.
BODY_SYSTEMS = [
    ("cardiovascular", {"en": "Cardiovascular", "es": "Cardiovascular", "icon": "\U0001FAC0"}),
    ("oncology", {"en": "Oncology / Cancer Risk", "es": "Oncología / Riesgo de Cáncer", "icon": "\U0001F397"}),
    ("neurological", {"en": "Neurological & Neurodevelopmental", "es": "Neurológico y Neurodesarrollo", "icon": "\U0001F9E0"}),
    ("metabolic_endocrine", {"en": "Metabolic & Endocrine", "es": "Metabólico y Endocrino", "icon": "⚖️"}),
    ("musculoskeletal", {"en": "Musculoskeletal & Connective Tissue", "es": "Musculoesquelético y Tejido Conectivo", "icon": "\U0001F9B4"}),
    ("hematologic_immune", {"en": "Hematologic & Immunological", "es": "Hematológico e Inmunológico", "icon": "\U0001FA78"}),
    ("respiratory", {"en": "Respiratory", "es": "Respiratorio", "icon": "\U0001FAC1"}),
    ("digestive_hepatic", {"en": "Digestive & Hepatic", "es": "Digestivo y Hepático", "icon": "\U0001F37D️"}),
    ("renal_urinary", {"en": "Renal & Urinary", "es": "Renal y Urinario", "icon": "\U0001FAD8"}),
    ("dermatologic_sensory", {"en": "Dermatologic, Eye & Ear", "es": "Dermatológico, Ojo y Oído", "icon": "\U0001F441️"}),
    ("reproductive_developmental", {"en": "Reproductive & Developmental", "es": "Reproductivo y del Desarrollo", "icon": "\U0001F476"}),
    ("other", {"en": "Other / Multisystem", "es": "Otro / Multisistémico", "icon": "\U0001F9EC"}),
]
BODY_SYSTEM_ORDER = [key for key, _ in BODY_SYSTEMS]
BODY_SYSTEM_LABELS = dict(BODY_SYSTEMS)

# Gene symbol (uppercase) -> system key. Deliberately partial: standard
# clinical-genetics gene-disease associations (OMIM/GeneCards-level
# knowledge), not individually citation-verified the way PRS panel PMIDs
# are - this is categorization, not a quantitative claim.
GENE_TO_SYSTEM = {
    # Cardiovascular
    "TBX5": "cardiovascular", "MYH7": "cardiovascular", "MYBPC3": "cardiovascular",
    "TNNT2": "cardiovascular", "TNNI3": "cardiovascular", "LMNA": "cardiovascular",
    "SCN5A": "cardiovascular", "KCNQ1": "cardiovascular", "KCNH2": "cardiovascular",
    "RYR2": "cardiovascular", "PKP2": "cardiovascular", "DSP": "cardiovascular",
    "DSG2": "cardiovascular", "FBN1": "cardiovascular", "COL3A1": "cardiovascular",
    "LDLR": "cardiovascular", "APOB": "cardiovascular", "PCSK9": "cardiovascular",
    "NOS3": "cardiovascular", "LRP8": "cardiovascular", "CNN2": "cardiovascular",
    "ECE1": "cardiovascular", "LGALS2": "cardiovascular", "CDKN2B": "cardiovascular",
    "CDKN2B-AS1": "cardiovascular", "F5": "cardiovascular", "F2": "cardiovascular",
    "MTHFR": "cardiovascular",
    # Oncology
    "BRCA1": "oncology", "BRCA2": "oncology", "TP53": "oncology", "RB1": "oncology",
    "VHL": "oncology", "RET": "oncology", "MEN1": "oncology", "NF1": "oncology",
    "NF2": "oncology", "MLH1": "oncology", "MSH2": "oncology", "MSH6": "oncology",
    "PMS2": "oncology", "APC": "oncology", "PTEN": "oncology", "STK11": "oncology",
    "CDH1": "oncology", "PALB2": "oncology", "ATM": "oncology", "CHEK2": "oncology",
    "LMO1": "oncology", "MYC": "oncology", "BMP4": "oncology",
    # Neurological / neurodevelopmental / psychiatric
    "IRF2BPL": "neurological", "ADAR": "neurological", "GRIN1": "neurological",
    "KCNB2": "neurological", "CNTNAP2": "neurological", "WDR45": "neurological",
    "CHAT": "neurological", "GJB1": "neurological", "CYP46A1": "neurological",
    "IBA57": "neurological", "APOE": "neurological", "PSEN1": "neurological",
    "PSEN2": "neurological", "APP": "neurological", "SNCA": "neurological",
    "LRRK2": "neurological", "MECP2": "neurological", "FMR1": "neurological",
    "SCN1A": "neurological", "COMT": "neurological", "BDNF": "neurological",
    "FOXO3": "neurological", "TOMM40": "neurological",
    "GABRA2": "neurological", "FKBP5": "neurological", "OPRM1": "neurological",
    "FAAH": "neurological", "ADRA2A": "metabolic_endocrine", "BMP2": "musculoskeletal",
    "CISH": "hematologic_immune", "MAPKAPK3": "hematologic_immune",
    "HS3ST1": "cardiovascular",
    # Metabolic / endocrine
    "NAGS": "metabolic_endocrine", "AQP7": "metabolic_endocrine", "CAPN10": "metabolic_endocrine",
    "TCF7L2": "metabolic_endocrine", "FTO": "metabolic_endocrine", "IGF2BP2": "metabolic_endocrine",
    "CDKAL1": "metabolic_endocrine", "VDR": "metabolic_endocrine", "CYP27B1": "metabolic_endocrine",
    "TG": "metabolic_endocrine", "ESR1": "metabolic_endocrine", "NAGLU": "metabolic_endocrine",
    "GAA": "metabolic_endocrine", "OTC": "metabolic_endocrine", "ATP7B": "metabolic_endocrine",
    "HFE": "metabolic_endocrine", "TTR": "metabolic_endocrine", "SOD2": "metabolic_endocrine",
    "PAH": "metabolic_endocrine", "G6PD": "metabolic_endocrine",
    # Musculoskeletal / connective tissue
    "RUNX2": "musculoskeletal", "CCN6": "musculoskeletal", "MMP3": "musculoskeletal",
    "CILP": "musculoskeletal", "COL1A1": "musculoskeletal", "COL5A1": "musculoskeletal",
    "TBXT": "musculoskeletal", "TFAP2A": "musculoskeletal", "RYR1": "musculoskeletal",
    "CACNA1S": "musculoskeletal", "ACTN3": "musculoskeletal",
    # Hematologic / immunological
    "IL10": "hematologic_immune", "IL19": "hematologic_immune", "IL13": "hematologic_immune",
    "IL1B": "hematologic_immune", "IL6": "hematologic_immune", "IL6-AS1": "hematologic_immune",
    "IRF4": "hematologic_immune", "IRF5": "hematologic_immune", "NLRP3": "hematologic_immune",
    "CCR2": "hematologic_immune", "CD209": "hematologic_immune", "CD244": "hematologic_immune",
    "SLC11A1": "hematologic_immune", "FCN3": "hematologic_immune", "MIF": "hematologic_immune",
    "MIF-AS1": "hematologic_immune", "MPO": "hematologic_immune", "EPO": "hematologic_immune",
    "GP1BB": "hematologic_immune", "ART4": "hematologic_immune", "CHI3L1": "hematologic_immune",
    "IL1RN": "hematologic_immune", "TNF": "hematologic_immune",
    # Respiratory
    "CFTR": "respiratory", "CFTR-AS1": "respiratory", "MUC5B": "respiratory",
    "NKX2-1": "respiratory", "SFTA3": "respiratory", "SCGB1A1": "respiratory",
    "HYKK": "respiratory", "PTGER2": "respiratory",
    # Digestive / hepatic
    "PRSS1": "digestive_hepatic", "TRB": "digestive_hepatic", "SPINK1": "digestive_hepatic",
    "HNF1A": "digestive_hepatic",
    # Renal / urinary
    "PKD1": "renal_urinary", "PKD2": "renal_urinary", "SLC9B1": "renal_urinary",
    "UMOD": "renal_urinary",
    # Dermatologic / eye / ear
    "KRT75": "dermatologic_sensory", "HERC2": "dermatologic_sensory", "OCA2": "dermatologic_sensory",
    "SLC45A2": "dermatologic_sensory", "SLC24A5": "dermatologic_sensory", "KITLG": "dermatologic_sensory",
    "TYR": "dermatologic_sensory", "MC1R": "dermatologic_sensory",
    # Reproductive / developmental
    "RNF212": "reproductive_developmental", "STOX1": "reproductive_developmental",
    "TUBB8": "reproductive_developmental",
}

# Substring (lowercase) -> system key, checked in order against the
# MedGen disease name when the gene lookup misses. Ordered so more
# specific terms are checked before generic ones (e.g. "renal cell
# carcinoma" should hit oncology before "renal").
_DISEASE_KEYWORDS = [
    ("carcinoma", "oncology"), ("neoplas", "oncology"), ("cancer", "oncology"),
    ("tumor", "oncology"), ("melanoma", "oncology"), ("leukemia", "oncology"),
    ("lymphoma", "oncology"), ("sarcoma", "oncology"),
    ("cardio", "cardiovascular"), ("heart", "cardiovascular"), ("aortic", "cardiovascular"),
    ("arrhythmia", "cardiovascular"), ("qt syndrome", "cardiovascular"),
    ("hypercholesterolemia", "cardiovascular"), ("thrombo", "cardiovascular"),
    ("neurodevelopmental", "neurological"), ("epilep", "neurological"), ("seizure", "neurological"),
    ("parkinson", "neurological"), ("alzheimer", "neurological"), ("dementia", "neurological"),
    ("charcot-marie", "neurological"), ("ataxia", "neurological"), ("neuropathy", "neurological"),
    ("intellectual disab", "neurological"), ("autism", "neurological"),
    ("diabetes", "metabolic_endocrine"), ("thyroid", "metabolic_endocrine"),
    ("hyperammonemia", "metabolic_endocrine"), ("metabolic", "metabolic_endocrine"),
    ("obesity", "metabolic_endocrine"), ("vitamin d", "metabolic_endocrine"),
    ("mucopolysaccharidosis", "metabolic_endocrine"), ("glycogen storage", "metabolic_endocrine"),
    ("dysplasia", "musculoskeletal"), ("osteogenesis", "musculoskeletal"),
    ("arthritis", "musculoskeletal"), ("myopathy", "musculoskeletal"), ("muscular dystrophy", "musculoskeletal"),
    ("holt-oram", "cardiovascular"), ("marfan", "cardiovascular"),
    ("cystic fibrosis", "respiratory"), ("pulmonary", "respiratory"), ("lung", "respiratory"),
    ("asthma", "respiratory"),
    ("pancreatitis", "digestive_hepatic"), ("crohn", "digestive_hepatic"), ("liver", "digestive_hepatic"),
    ("hepatic", "digestive_hepatic"), ("colitis", "digestive_hepatic"),
    ("kidney", "renal_urinary"), ("renal", "renal_urinary"), ("nephro", "renal_urinary"),
    ("usher syndrome", "dermatologic_sensory"), ("retinitis", "dermatologic_sensory"),
    ("hearing loss", "dermatologic_sensory"), ("deafness", "dermatologic_sensory"),
    ("blindness", "dermatologic_sensory"), ("skin", "dermatologic_sensory"),
    ("albinism", "dermatologic_sensory"), ("ichthyosis", "dermatologic_sensory"),
    ("hemophilia", "hematologic_immune"), ("anemia", "hematologic_immune"),
    ("immunodeficiency", "hematologic_immune"), ("leprosy", "hematologic_immune"),
    ("hepatitis", "hematologic_immune"),
    ("preeclampsia", "reproductive_developmental"), ("infertility", "reproductive_developmental"),
    ("syndrome", "other"),  # generic eponymous syndromes with no other hit - last resort before "other"
]


def _clean_gene_symbol(raw: str) -> str:
    """gene_info fields look like 'CFTR:1080|CFTR-AS1:111082987' - take the first symbol."""
    if not raw:
        return ""
    first = raw.split("|")[0].split(":")[0].strip().upper()
    return first


def classify_body_system(genes, disease_name: str = "") -> str:
    """
    Classify a ClinVar variant into a body system.

    `genes` may be a list of gene symbols (preferred, e.g. the `genes`
    field already parsed by clinvar_annotator.py) or a raw `gene_info`
    string like "CFTR:1080|CFTR-AS1:111082987". Falls back to keyword
    matching on `disease_name`, then to "other".
    """
    if isinstance(genes, str):
        genes = [_clean_gene_symbol(genes)]
    for g in genes or []:
        symbol = _clean_gene_symbol(g) if ":" in (g or "") or "|" in (g or "") else (g or "").upper()
        if symbol in GENE_TO_SYSTEM:
            return GENE_TO_SYSTEM[symbol]

    name = (disease_name or "").lower().replace("_", " ")
    if name and name not in (".", "not provided"):
        for keyword, system in _DISEASE_KEYWORDS:
            if keyword in name:
                return system

    return "other"


def system_label(system_key: str, lang: str = "en") -> str:
    entry = BODY_SYSTEM_LABELS.get(system_key, BODY_SYSTEM_LABELS["other"])
    return f"{entry['icon']} {entry.get(lang, entry['en'])}"
