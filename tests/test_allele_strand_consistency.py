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


# ── Palindromic (A/T or C/G) SNP strand-orientation tracking ──────────────
#
# A/T and C/G allele pairs are self-complementary: the panel's stated
# effect_allele/reference_allele pair looks identical whether the source
# data was truly forward-strand or was silently flipped. The checks above
# are structurally blind to this - a flipped palindromic SNP still
# "matches" the fixture's allele set. There is no risk_allele_frequency/
# MAF column in snp_database_annotated.csv for a generic quantitative
# cross-check against 1000G population frequencies (documented
# limitation, out of scope here). This section does not attempt to
# resolve the unresolved rows below - it makes the existing risk visible
# and gates future regressions (IMPROVEMENT_PLAN.md 4.3/4.1).

PALINDROMIC_PAIRS = {frozenset({"A", "T"}), frozenset({"C", "G"})}


def _is_palindromic(effect_allele, reference_allele):
    ea, ra = effect_allele.strip(), reference_allele.strip()
    return len(ea) == 1 and len(ra) == 1 and frozenset({ea, ra}) in PALINDROMIC_PAIRS


# Palindromic rows whose forward-strand orientation IS independently
# documented in the CSV's notes column (dbSNP gene-strand annotation,
# Ensembl VEP coding strand, or a 1000G population-frequency
# cross-check) - not just letter-membership.
PALINDROME_STRAND_RESOLVED = {
    ("rs1799945", "Iron levels"):
        "notes: verified against dbSNP GRCh37 (chr6:26091179, C>G); HFE is "
        "gene-plus-strand, so gene-relative and genomic-forward notation "
        "agree here (2026-07-14).",
    ("rs16891982", "Hair color (black)"):
        "notes: verified directly against local 1000G reference genotypes "
        "(G=93.8% EUR vs 3.6% AFR) - G is the derived light-pigmentation "
        "allele, matching its use elsewhere in this panel (2026-08-06).",
    ("rs16891982", "Eye color"):
        "notes: coding direction verified via Ensembl VEP (minus strand) "
        "+ the same 1000G population-frequency cross-check as the Hair "
        "color (black) row above (2026-08-06).",
    ("rs16891982", "Skin pigmentation"):
        "same physical SNP/position as the Hair color (black)/Eye color "
        "rows above (effect_allele=G/reference_allele=C here, consistent "
        "with that resolution).",
    ("rs16891982", "Hair color (brown)"):
        "same physical SNP/position as the Hair color (black)/Eye color "
        "rows above (effect_allele=C/reference_allele=G here, the "
        "complementary direction for the dark-hair effect allele, "
        "consistent with that resolution).",
}

# Palindromic rows with NO independent strand-orientation verification on
# record - a real, currently-untracked risk. To clear one:
#   1. Verify its true forward-strand orientation against an independent
#      source (dbSNP strand annotation, Ensembl VEP, or a population
#      allele frequency in reference/1000G_full matching a literature-
#      reported frequency for a named allele) - letter-matching against
#      the panel's own data is NOT sufficient for a palindromic pair.
#   2. Record a dated note in snp_database_annotated.csv's `notes` column
#      describing how it was verified.
#   3. Move the entry here to PALINDROME_STRAND_RESOLVED citing that note.
# This list should shrink over time, not grow.
#
# (rs6968865's note only confirms a corrected letter is *a member* of the
# bim's observed pair - exactly the palindrome-blind check, not an
# independent resolution. rs2304672's note fixes which allele is the
# literature-reported morning-preference allele - a phenotype-direction
# fix, not a forward-strand cross-check. Both deliberately stay here.)
KNOWN_UNRESOLVED_PALINDROMES = {
    ("rs6968865", "Caffeine metabolism"),
    ("rs328", "Lipid metabolism"),
    ("rs3135506", "Lipid metabolism"),
    ("rs964184", "Lipid metabolism"),
    ("rs174548", "Omega-3 metabolism"),
    ("rs174575", "Omega-3 metabolism"),
    ("rs2236212", "Omega-3 metabolism"),
    ("rs145946881", "Lactose intolerance"),
    ("rs41525747", "Lactose intolerance"),
    ("rs9939609", "Obesity predisposition"),
    ("rs1558902", "Obesity predisposition"),
    ("rs10767664", "Obesity predisposition"),
    ("rs4818", "Dopamine regulation"),
    ("rs1042522", "Detoxification"),
    ("rs2910164", "Detoxification"),
    ("rs6013897", "Vitamin D metabolism"),
    ("rs4506565", "Glucose metabolism"),
    ("rs10830963", "Glucose metabolism"),
    ("rs1801253", "Blood pressure"),
    ("rs713598", "Bitter taste perception"),
    ("rs9939609", "Body composition"),
    ("rs2304672", "Morning chronotype (early bird)"),
}


def test_palindromic_snps_are_explicitly_triaged():
    """Every A/T or C/G effect/reference-allele pair in the panel must be
    in PALINDROME_STRAND_RESOLVED (independently verified) or
    KNOWN_UNRESOLVED_PALINDROMES (tracked gap). A new palindromic SNP
    landing untriaged fails CI by design, forcing the same manual
    verification already applied to other rows in this panel's history
    instead of silently inheriting an invisible strand-flip risk."""
    rows = _load_panel_rows()
    untracked = []
    for r in rows:
        if not _is_palindromic(r["effect_allele"], r["reference_allele"]):
            continue
        key = (r["rsid"], r["trait_category"])
        if key in PALINDROME_STRAND_RESOLVED or key in KNOWN_UNRESOLVED_PALINDROMES:
            continue
        untracked.append(
            f"{r['rsid']} ({r['trait_category']}): "
            f"{r['effect_allele']}/{r['reference_allele']} is palindromic "
            f"and not triaged"
        )
    assert not untracked, (
        "New palindromic SNP(s) with no strand-orientation triage - A/T "
        "and C/G pairs are self-complementary, so the usual forward-strand "
        "letter check can't tell whether they were read from the correct "
        "genomic strand:\n" + "\n".join(untracked) +
        "\n\nTo fix: verify the true forward-strand orientation (dbSNP, "
        "Ensembl VEP, or a 1000G population allele frequency matching a "
        "literature-reported frequency), note it in snp_database_annotated"
        ".csv, and add to PALINDROME_STRAND_RESOLVED citing that note - or "
        "add to KNOWN_UNRESOLVED_PALINDROMES to track the gap explicitly."
    )


def test_palindrome_allow_lists_have_no_stale_entries():
    """Guards against drift: every listed (rsid, trait_category) must
    still be an actual palindromic row in the current panel."""
    rows = _load_panel_rows()
    current = {
        (r["rsid"], r["trait_category"]) for r in rows
        if _is_palindromic(r["effect_allele"], r["reference_allele"])
    }
    listed = set(PALINDROME_STRAND_RESOLVED) | KNOWN_UNRESOLVED_PALINDROMES
    stale = listed - current
    assert not stale, f"Stale palindrome allow-list entries: {stale}"
