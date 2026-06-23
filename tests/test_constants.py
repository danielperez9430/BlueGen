"""Tests for utils.constants module."""

import sys
from pathlib import Path

# Add scripts to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from utils.constants import (
    PIPELINE_DIR,
    PROJECT_ROOT,
    SCRIPTS_DIR,
    TOOLS_DIR,
    FILE_SNP_DB,
    TRAIT_CATEGORIES,
    RISK_LOW,
    RISK_HIGH,
    PLINK_THREADS,
    PLINK_MEMORY,
    PIPELINE_VERSION,
    GENOME_BUILD,
)


def test_pipeline_dir_exists():
    assert PIPELINE_DIR.exists()
    assert PIPELINE_DIR.is_dir()


def test_project_root_is_parent():
    assert PROJECT_ROOT == PIPELINE_DIR.parent


def test_scripts_dir_exists():
    assert SCRIPTS_DIR.exists()
    assert SCRIPTS_DIR.is_dir()


def test_tools_dir_path():
    assert "tools" in str(TOOLS_DIR)


def test_snp_db_is_string():
    assert isinstance(FILE_SNP_DB, str)
    assert FILE_SNP_DB.endswith(".csv")


def test_trait_categories():
    assert isinstance(TRAIT_CATEGORIES, list)
    assert len(TRAIT_CATEGORIES) > 0
    assert all(isinstance(t, str) for t in TRAIT_CATEGORIES)


def test_risk_thresholds():
    assert RISK_LOW < RISK_HIGH
    assert 0 <= RISK_LOW <= 50
    assert 50 <= RISK_HIGH <= 100


def test_plink_defaults():
    assert PLINK_THREADS > 0
    assert PLINK_MEMORY > 0


def test_pipeline_version():
    assert isinstance(PIPELINE_VERSION, str)
    assert "." in PIPELINE_VERSION


def test_genome_build():
    assert GENOME_BUILD in ("GRCh37", "GRCh38")
