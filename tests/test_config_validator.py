"""Tests for utils.config_validator module."""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from utils.config_validator import validate_config, SCHEMA


def test_validate_actual_config():
    """Validate the real config.yaml."""
    config_path = Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "config.yaml"
    if config_path.exists():
        config = validate_config(str(config_path))
        assert isinstance(config, dict)
        assert "pipeline" in config


def test_validate_missing_file_exits():
    """Missing config file should raise SystemExit."""
    import pytest
    with pytest.raises(SystemExit):
        validate_config("/nonexistent/config.yaml")


def test_validate_invalid_yaml_exits():
    """Invalid YAML should raise SystemExit."""
    import pytest
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("{{{{invalid yaml")
        f.flush()
        try:
            with pytest.raises(SystemExit):
                validate_config(f.name)
        finally:
            Path(f.name).unlink()


def test_validate_empty_config_exits():
    """Empty config should raise SystemExit (missing required fields)."""
    import pytest
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("pipeline:\n  name: test\n")
        f.flush()
        try:
            with pytest.raises(SystemExit):
                validate_config(f.name)
        finally:
            Path(f.name).unlink()


def test_schema_has_required_sections():
    """Schema should define all major pipeline sections."""
    required_sections = ["pipeline", "input", "plink", "qc", "prs", "logging"]
    for section in required_sections:
        assert section in SCHEMA, f"Missing schema section: {section}"


def test_validate_good_config():
    """A well-formed config should pass validation."""
    import pytest
    import yaml
    good_config = {
        "pipeline": {"name": "test", "version": "1.0", "genome_build": "GRCh37"},
        "input": {"vcf": "test.vcf.gz", "sample_id": "TEST"},
        "plink": {"threads": 4, "memory": 8000},
        "qc": {"snp_missingness": 0.1, "individual_missingness": 0.1, "maf": 0.01, "hwe": 0.000001},
        "ld_pruning": {"r2_threshold": 0.2},
        "pca": {"num_components": 10},
        "prs": {"trait_categories": ["test"], "score_method": "sum", "risk_thresholds": {"low": 25, "high": 75}},
        "population_calibration": {"enabled": True, "num_pcs": 20},
        "normalization": {"method": "percentile"},
        "bilingual": {"enabled": True},
        "logging": {"level": "INFO"},
    }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(good_config, f)
        f.flush()
        try:
            config = validate_config(f.name)
            assert config["pipeline"]["name"] == "test"
        finally:
            Path(f.name).unlink()
