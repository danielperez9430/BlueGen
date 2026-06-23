#!/usr/bin/env python3
"""
Centralized constants for the PRS Research Pipeline.

Eliminates magic strings and duplicated path definitions across
prs.py, run.sh, and individual scripts.

Usage:
    from utils.constants import PIPELINE_DIR, TRAIT_CATEGORIES, SNP_DB_PATH
"""

from pathlib import Path

# ── Project paths ─────────────────────────────────────────────────────────────

PIPELINE_DIR = Path(__file__).resolve().parent.parent.parent  # utils/ → scripts/ → prs_research_pipeline/
PROJECT_ROOT = PIPELINE_DIR.parent
SCRIPTS_DIR = PIPELINE_DIR / "scripts"
TOOLS_DIR = PROJECT_ROOT / "tools"

# ── Output directories (relative to PIPELINE_DIR) ────────────────────────────

DIR_PLINK = "plink"
DIR_QC = "qc"
DIR_PCA = "pca"
DIR_PRS = "prs"
DIR_SCIENCE = "science"
DIR_BENCHMARK = "benchmark"
DIR_REPORTS = "reports"
DIR_CLINVAR = "clinvar"
DIR_PHARMGKB = "pharmgkb"
DIR_ANCESTRY = "ancestry"
DIR_INTERPRETATIONS = "interpretations"
DIR_REPRODUCIBILITY = "reproducibility"
DIR_VALIDATION = "validation"
DIR_PUBLICATION = "publication_evidence_pack"

# ── Key output files (relative to PIPELINE_DIR) ──────────────────────────────

FILE_SNP_DB = "data/snp_database_annotated.csv"
FILE_PRS_RAW = "prs/prs_raw.csv"
FILE_PRS_ADJUSTED = "prs/prs_adjusted.csv"
FILE_PRS_CALIBRATED_V2 = "prs/population_calibrated_v2.csv"
FILE_PRS_CALIBRATED = "prs/population_calibrated.csv"
FILE_PRS_RESULT = "prs/PRS_RESULT.json"
FILE_ANCESTRY_MODEL = "science/ANCESTRY_MODEL.json"
FILE_ANCESTRY_CLASSIFICATION = "ancestry/classification_report.json"
FILE_ANCESTRY_LEGACY = "pca/ancestry_inference.json"
FILE_TARGET_PCS = "pca/target_pcs.eigenvec"
FILE_1000G_PCS = "pca/1000G_pca.eigenvec"
FILE_VALIDATION_REPORT = "science/global_validation_report.json"
FILE_FINAL_SCORE = "FINAL_SCIENTIFIC_SCORE.json"
FILE_PUBLICATION_LOCK = "PUBLICATION_LOCK.md"

# ── Reference data paths ─────────────────────────────────────────────────────

REF_1000G_FULL_DIR = "reference/1000G_full"
REF_1000G_FULL_BFILE = "reference/1000G_full/1000G_full"
REF_1000G_CHR22_VCF = "reference/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
REF_POP_PANEL_FULL = "reference/1000G_full/population_panel.txt"
REF_POP_PANEL_CHR22 = "reference/1000G/20130606_g1k_3202_samples_ped_population.txt"
REF_CLINVAR_VCF = "reference/clinvar/clinvar.vcf.gz"
REF_MEDGEN_DIR = "reference/medgen"
REF_CLINPGX_DIR = "reference/clinpgx"

# ── PRS trait categories ─────────────────────────────────────────────────────

TRAIT_CATEGORIES = [
    "Caffeine metabolism",
    "Lipid metabolism",
    "Omega-3 metabolism",
    "Folate & methylation",
    "Lactose intolerance",
    "Obesity predisposition",
    "Dopamine regulation",
    "Detoxification",
    "Vitamin D metabolism",
    "Glucose metabolism",
]

# ── Risk thresholds ──────────────────────────────────────────────────────────

RISK_LOW = 25
RISK_HIGH = 75

# ── PLINK defaults ───────────────────────────────────────────────────────────

PLINK_THREADS = 4
PLINK_MEMORY = 8000

# ── Version ──────────────────────────────────────────────────────────────────

PIPELINE_VERSION = "2.0.0"
GENOME_BUILD = "GRCh37"
