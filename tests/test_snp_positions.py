"""Regression tests for the curated SNP panel's GRCh37 positions.

Wraps scripts/setup/audit_snp_positions.py (plan tasks 3.3 / H4): the PRS
scorer matches SNPs by chrom:pos (prs_plink_score.py), so a wrong position
makes a SNP silently fail to match in PLINK and its trait drops out of the
report with no error. These tests catch that class of bug before it reaches
main.
"""

import csv
import sys
import urllib.error
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from setup.audit_snp_positions import DEFAULT_CSV, ensembl_batch, grch37_mapping, load_rows, norm_chrom

CANONICAL_CHROMS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}


def _rows():
    _, rows = load_rows(DEFAULT_CSV)
    return rows


def test_csv_exists_and_loads():
    assert DEFAULT_CSV.exists()
    rows = _rows()
    assert len(rows) > 0


def test_every_row_has_a_canonical_chrom_and_positive_pos():
    for r in _rows():
        rsid = r["rsid"].strip()
        chrom = norm_chrom(r.get("chrom", ""))
        assert chrom in CANONICAL_CHROMS, f"{rsid}: non-canonical chrom {r.get('chrom')!r}"
        pos = r.get("pos", "").strip()
        assert pos.isdigit() and int(pos) > 0, f"{rsid}: invalid pos {pos!r}"


def test_every_row_has_rsid_and_alleles():
    for r in _rows():
        rsid = r["rsid"].strip()
        assert rsid.startswith("rs"), f"malformed rsid: {rsid!r}"
        assert r.get("effect_allele", "").strip(), f"{rsid}: empty effect_allele"
        assert r.get("reference_allele", "").strip(), f"{rsid}: empty reference_allele"


def test_no_position_drift_against_ensembl_grch37():
    """The real regression check: every rsID's CSV position must still match
    Ensembl GRCh37. Skips (rather than fails) if Ensembl is unreachable, so a
    transient network/API outage doesn't block unrelated merges - but any
    actual chrom/pos mismatch fails the build.
    """
    rows = _rows()
    rsids = sorted({r["rsid"].strip() for r in rows if r.get("rsid", "").strip()})

    try:
        data = ensembl_batch(rsids)
    except (urllib.error.URLError, RuntimeError, TimeoutError) as e:
        pytest.skip(f"Ensembl GRCh37 REST unreachable, skipping live position check: {e}")

    mismatches = []
    for r in rows:
        rsid = r["rsid"].strip()
        v = data.get(rsid)
        if not v:
            mismatches.append(f"{rsid}: not found in Ensembl GRCh37")
            continue
        ref = grch37_mapping(v)
        if ref is None:
            mismatches.append(f"{rsid}: no canonical GRCh37 mapping returned")
            continue
        ens_chrom, ens_pos = ref
        csv_chrom, csv_pos = norm_chrom(r.get("chrom", "")), r.get("pos", "").strip()
        if csv_chrom != ens_chrom:
            mismatches.append(
                f"{rsid}: CHROM ERROR csv=chr{csv_chrom}:{csv_pos} ensembl=chr{ens_chrom}:{ens_pos}"
            )
        elif csv_pos != str(ens_pos):
            mismatches.append(
                f"{rsid}: POS ERROR csv=chr{csv_chrom}:{csv_pos} ensembl=chr{ens_chrom}:{ens_pos}"
            )

    assert not mismatches, (
        f"{len(mismatches)} SNP(s) diverge from Ensembl GRCh37 - these will silently "
        f"fail to match in PLINK --score:\n" + "\n".join(mismatches)
    )


def test_no_duplicate_rsid_trait_pairs():
    """Catches accidental copy-paste duplicate rows for the same trait."""
    rows = _rows()
    seen = {}
    dupes = []
    for r in rows:
        key = (r["rsid"].strip(), r["trait_category"].strip())
        if key in seen:
            dupes.append(key)
        seen[key] = True
    assert not dupes, f"duplicate (rsid, trait_category) rows: {dupes}"


def test_csv_matches_audit_script_default_path():
    """Guards against the CSV moving without audit_snp_positions.py being updated."""
    with open(DEFAULT_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
    for col in ("rsid", "gene", "trait_category", "effect_allele", "reference_allele", "chrom", "pos"):
        assert col in fieldnames, f"expected column {col!r} missing from {DEFAULT_CSV.name}"
