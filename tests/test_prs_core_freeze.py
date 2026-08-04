"""Unit tests for PRS_CORE freeze/refresh logic (scripts/sss/36_prs_core.py).

Regression coverage for the bug where prs_core_definition.json was frozen
once (2026-06-03, 109 SNPs/10 traits) and never refreshed despite the SNP
panel being curated up to 190 SNPs/59 traits over the following weeks -
silently leaving the canonical SSST definition, and everything that reads
it (manuscripts, consistency checks), describing a stale panel.
"""

import csv
import importlib.util
import os

import pytest

_MODULE_PATH = os.path.join(os.path.dirname(__file__), "..",
    "prs_research_pipeline", "scripts", "sss", "36_prs_core.py")
_spec = importlib.util.spec_from_file_location("prs_core_module", _MODULE_PATH)
_prs_core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_prs_core)
PRSCoreRegistry = _prs_core.PRSCoreRegistry


def _write_csv(path, rows):
    fieldnames = ["rsid", "gene", "trait_category", "effect_allele",
                  "reference_allele", "weight", "evidence_level", "notes",
                  "pmid", "chrom", "pos"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _row(rsid, trait, gene="GENE"):
    return {"rsid": rsid, "gene": gene, "trait_category": trait,
            "effect_allele": "A", "reference_allele": "G", "weight": "0.2",
            "evidence_level": "B", "notes": "", "pmid": "12345",
            "chrom": "chr1", "pos": "1000"}


def test_creates_fresh_definition_with_csv_hash(tmp_path):
    snp_db = tmp_path / "panel.csv"
    _write_csv(snp_db, [_row("rs1", "Trait A"), _row("rs2", "Trait B")])

    registry = PRSCoreRegistry(output_dir=str(tmp_path / "science"))
    core = registry.load_or_create(snp_db=str(snp_db))

    assert core.n_variants == 2
    assert core.n_traits == 2
    assert core.source_csv_hash != ""


def test_unchanged_csv_reuses_frozen_definition(tmp_path):
    snp_db = tmp_path / "panel.csv"
    _write_csv(snp_db, [_row("rs1", "Trait A")])

    registry = PRSCoreRegistry(output_dir=str(tmp_path / "science"))
    first = registry.load_or_create(snp_db=str(snp_db))

    # Second call, same CSV content, same registry instance: must not re-freeze.
    second = registry.load_or_create(snp_db=str(snp_db))
    assert second.frozen_date == first.frozen_date
    assert second.source_csv_hash == first.source_csv_hash


def test_changed_csv_triggers_refresh(tmp_path):
    snp_db = tmp_path / "panel.csv"
    _write_csv(snp_db, [_row("rs1", "Trait A")])

    registry = PRSCoreRegistry(output_dir=str(tmp_path / "science"))
    stale = registry.load_or_create(snp_db=str(snp_db))
    assert stale.n_variants == 1

    # Panel grows - this is what H1-H4-style curation does to the real CSV.
    _write_csv(snp_db, [_row("rs1", "Trait A"), _row("rs2", "Trait B"),
                         _row("rs3", "Trait C")])

    refreshed = registry.load_or_create(snp_db=str(snp_db))
    assert refreshed.n_variants == 3
    assert refreshed.n_traits == 3
    assert refreshed.source_csv_hash != stale.source_csv_hash


def test_legacy_frozen_file_without_hash_field_triggers_refresh(tmp_path):
    """Files frozen before this fix have no source_csv_hash key at all -
    _load() defaults it to "", which must not equal a real hash and so
    must not be treated as 'unchanged'."""
    snp_db = tmp_path / "panel.csv"
    _write_csv(snp_db, [_row("rs1", "Trait A"), _row("rs2", "Trait B")])

    science_dir = tmp_path / "science"
    science_dir.mkdir()
    legacy_json = science_dir / "prs_core_definition.json"
    import json
    legacy_json.write_text(json.dumps({
        "formula": "PRS = Σ(βⱼ × Gᵢⱼ)", "n_variants": 109, "n_traits": 10,
        "traits": ["Old Trait"], "frozen_date": "2026-06-03 12:45 UTC",
        "definition_hash": "deadbeef",
        # no source_csv_hash key - simulates a pre-fix frozen file
    }))

    registry = PRSCoreRegistry(output_dir=str(science_dir))
    core = registry.load_or_create(snp_db=str(snp_db))

    assert core.n_variants == 2
    assert core.n_traits == 2
