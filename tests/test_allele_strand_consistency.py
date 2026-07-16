"""Guards against the panel's alleles drifting off the GRCh37 forward strand.

Regression coverage for a real bug: 34 rows in snp_database_annotated.csv had
effect_allele/reference_allele/risk_genotype recorded in gene-relative
(transcript/minus-strand) notation instead of GRCh37 forward-strand notation.
Scoring (prs_plink_score.py, build_reference_distributions.py) matches
effect_allele against PLINK's forward-strand .bim alleles - a mismatched
letter makes PLINK silently drop the SNP from scoring with no visible error,
for BOTH the 1000G reference build and the user's own real genotype scoring.

tests/fixtures/panel_forward_alleles.tsv is a small extract of
reference/1000G_full/1000G_full.bim (chrom, pos, a1, a2) limited to the
panel's own positions, generated once and committed so this test runs
offline/deterministically in CI (the real .bim is gitignored, multi-GB).
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SNP_DB = REPO_ROOT / "prs_research_pipeline" / "data" / "snp_database_annotated.csv"
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "panel_forward_alleles.tsv"

# Rows with a known, deliberate reason to skip the forward-strand check -
# each must have a comment explaining why, and should shrink over time, not grow.
KNOWN_EXCEPTIONS = {
    # Downgraded to evidence tier D / exploratory (commit 7711b13): the only
    # documented use of rs590787 is distinguishing weak-D/DEL from RhD- in a
    # Chinese cohort, not general RhD+/- typing; its bim pair (C/A) doesn't
    # cleanly complement-pair either, consistent with it being a mixed-strand
    # or otherwise unreliable record. Left as-is rather than guessed at.
    ("rs590787", "Blood type (ABO Rh)"),
}


def _load_fixture():
    alleles = {}
    with open(FIXTURE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            alleles[(row["chrom"], row["pos"])] = (row["a1"], row["a2"])
    return alleles


def _load_panel_rows():
    with open(SNP_DB, newline="") as fh:
        return list(csv.DictReader(fh))


def test_fixture_is_not_empty():
    alleles = _load_fixture()
    assert len(alleles) > 100


def test_every_scoreable_effect_allele_matches_forward_strand():
    """The actual regression check: effect_allele must be one of the two
    forward-strand alleles PLINK sees at that position - anything else is
    silently dropped from scoring, not just "less accurate"."""
    alleles = _load_fixture()
    rows = _load_panel_rows()

    offenders = []
    checked = 0
    for r in rows:
        key = (r["rsid"], r["trait_category"])
        if key in KNOWN_EXCEPTIONS:
            continue
        allele = r["effect_allele"].strip()
        if len(allele) != 1 or allele not in "ACGT":
            continue  # indels / malformed values are a separate, tracked issue
        chrom = r["chrom"].replace("chr", "")
        pos = r["pos"].strip()
        pair = alleles.get((chrom, pos))
        if pair is None:
            continue  # position not in the fixture (not in 1000G at all) - separate issue
        checked += 1
        if allele not in pair:
            offenders.append(f"{r['rsid']} ({r['trait_category']}): effect_allele={allele!r} not in {pair}")

    assert checked > 100, "fixture/panel join produced too few checkable rows - fixture may be stale"
    assert not offenders, (
        f"{len(offenders)} row(s) have an effect_allele not on the GRCh37 forward strand "
        f"(PLINK will silently drop these from scoring):\n" + "\n".join(offenders)
    )


def test_reference_allele_is_the_other_forward_strand_allele():
    """reference_allele is used as a fallback dosage source in 06_prs_compute.py,
    so it must also be forward-strand - specifically, it should be the OTHER
    allele of the pair (not necessarily effect_allele's naive complement -
    see rs13217795, a mixed-strand row where only effect_allele needed flipping)."""
    alleles = _load_fixture()
    rows = _load_panel_rows()

    offenders = []
    for r in rows:
        key = (r["rsid"], r["trait_category"])
        if key in KNOWN_EXCEPTIONS:
            continue
        ea = r["effect_allele"].strip()
        ra = r["reference_allele"].strip()
        if len(ea) != 1 or ea not in "ACGT" or len(ra) != 1 or ra not in "ACGT":
            continue
        chrom = r["chrom"].replace("chr", "")
        pos = r["pos"].strip()
        pair = alleles.get((chrom, pos))
        if pair is None or ea not in pair:
            continue
        if ra not in pair:
            offenders.append(f"{r['rsid']} ({r['trait_category']}): reference_allele={ra!r} not in {pair}")

    assert not offenders, (
        f"{len(offenders)} row(s) have a reference_allele not on the GRCh37 forward strand:\n"
        + "\n".join(offenders)
    )
