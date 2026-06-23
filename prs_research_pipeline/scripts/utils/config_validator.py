#!/usr/bin/env python3
"""
Configuration schema validation for the PRS Research Pipeline.

Validates config.yaml at startup to fail fast on misconfiguration.
Checks types, ranges, and allowed values for all pipeline parameters.

Schema covers:
    - pipeline: name, version, genome_build
    - input: vcf, sample_id
    - plink: threads (1-64), memory (512-131072 MB)
    - qc: missingness thresholds, MAF, HWE
    - ld_pruning: r2_threshold (0-1)
    - prs: trait_categories, score_method, risk_thresholds
    - population_calibration: enabled, num_pcs
    - normalization: method (zscore|percentile|internal)
    - logging: level (DEBUG|INFO|WARNING|ERROR)

Usage:
    from utils.config_validator import validate_config
    config = validate_config("config.yaml")  # raises SystemExit on error

    # CLI
    python config_validator.py config.yaml
"""

import sys
from pathlib import Path
from typing import Any

import yaml


# ── Schema definition ─────────────────────────────────────────────────────────
# Each key maps to a validator: (required, type, [min, max] or [allowed_values])

SCHEMA = {
    "pipeline": {
        "name": (True, str),
        "version": (True, str),
        "genome_build": (True, str, ["GRCh37", "GRCh38"]),
    },
    "input": {
        "vcf": (True, str),
        "sample_id": (True, str),
    },
    "plink": {
        "threads": (True, int, 1, 64),
        "memory": (True, int, 512, 131072),
    },
    "qc": {
        "snp_missingness": (True, float, 0.0, 1.0),
        "individual_missingness": (True, float, 0.0, 1.0),
        "maf": (True, float, 0.0, 0.5),
        "hwe": (True, float, 0.0, 1.0),
    },
    "ld_pruning": {
        "r2_threshold": (True, float, 0.0, 1.0),
    },
    "pca": {
        "num_components": (True, int, 1, 100),
    },
    "prs": {
        "trait_categories": (True, list),
        "score_method": (True, str, ["sum", "avg", "std"]),
        "risk_thresholds": {
            "low": (True, int, 0, 50),
            "high": (True, int, 50, 100),
        },
    },
    "population_calibration": {
        "enabled": (True, bool),
        "num_pcs": (True, int, 1, 100),
    },
    "normalization": {
        "method": (True, str, ["zscore", "percentile", "internal"]),
    },
    "bilingual": {
        "enabled": (True, bool),
    },
    "logging": {
        "level": (True, str, ["DEBUG", "INFO", "WARNING", "ERROR"]),
    },
}


def _validate_value(value: Any, spec: tuple, path: str) -> list[str]:
    """Validate a single value against its spec. Returns list of error messages."""
    errors = []
    required = spec[0]
    expected_type = spec[1]

    if value is None:
        if required:
            errors.append(f"{path}: required field is missing")
        return errors

    if not isinstance(value, expected_type):
        errors.append(f"{path}: expected {expected_type.__name__}, got {type(value).__name__}")
        return errors

    # Range check (for int/float with min/max args)
    if len(spec) == 4 and expected_type in (int, float):
        min_val, max_val = spec[2], spec[3]
        if value < min_val or value > max_val:
            errors.append(f"{path}: value {value} out of range [{min_val}, {max_val}]")

    # Allowed values check
    if len(spec) == 3 and isinstance(spec[2], list):
        if value not in spec[2]:
            errors.append(f"{path}: value '{value}' not in allowed: {spec[2]}")

    return errors


def _validate_section(config: dict, schema: dict, prefix: str = "") -> list[str]:
    """Recursively validate a config section against schema."""
    errors = []

    for key, spec in schema.items():
        path = f"{prefix}.{key}" if prefix else key

        if isinstance(spec, dict):
            # Nested section
            value = config.get(key)
            if value is None:
                errors.append(f"{path}: required section is missing")
            elif not isinstance(value, dict):
                errors.append(f"{path}: expected section (dict), got {type(value).__name__}")
            else:
                errors.extend(_validate_section(value, spec, path))
        else:
            # Leaf value
            value = config.get(key)
            errors.extend(_validate_value(value, spec, path))

    return errors


def validate_config(config_path: str = "config.yaml") -> dict:
    """
    Validate config.yaml against the schema.

    Returns parsed config dict on success.
    Prints errors and calls sys.exit(1) on failure.
    """
    path = Path(config_path)
    if not path.exists():
        print(f"\033[0;31m✗ Config not found: {config_path}\033[0m")
        sys.exit(1)

    try:
        with open(path) as fh:
            config = yaml.safe_load(fh)
    except yaml.YAMLError as e:
        print(f"\033[0;31m✗ Invalid YAML: {e}\033[0m")
        sys.exit(1)

    if not isinstance(config, dict):
        print(f"\033[0;31m✗ Config must be a YAML mapping, got {type(config).__name__}\033[0m")
        sys.exit(1)

    errors = _validate_section(config, SCHEMA)

    if errors:
        print(f"\n\033[0;31m✗ Config validation failed ({len(errors)} errors):\033[0m")
        for e in errors:
            print(f"  • {e}")
        print(f"\n  Fix: {config_path}")
        sys.exit(1)

    return config


def main():
    """CLI entry point for config validation."""
    import argparse
    parser = argparse.ArgumentParser(description="Validate pipeline config.yaml")
    parser.add_argument("config", nargs="?", default="config.yaml", help="Config file path")
    args = parser.parse_args()

    config = validate_config(args.config)
    print(f"\033[0;32m✓ Config valid: {args.config}\033[0m")
    print(f"  Pipeline: {config.get('pipeline', {}).get('name', '?')}")
    print(f"  Version:  {config.get('pipeline', {}).get('version', '?')}")
    print(f"  Build:    {config.get('pipeline', {}).get('genome_build', '?')}")


if __name__ == "__main__":
    main()
