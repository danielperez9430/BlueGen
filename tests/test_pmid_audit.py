"""Literature-citation audit of the curated panel and recommendations
(RELEASE_PLAN 3.0.6), run OFFLINE against the committed PubMed cache.

Every PMID cited by data/snp_database_annotated.csv or
data/trait_recommendations.json must (a) be in tests/fixtures/pubmed_cache.json
- otherwise run `scripts/setup/audit_pmids.py --update-cache` and commit the
cache - and (b) have a PubMed title/abstract that mentions the rsID, a gene
of the row (or its protein name), or the trait, unless a reviewed reason is
recorded in data/pmid_audit_allowlist.json. On 2026-09-15 this check found
39 panel rows citing papers about mucormycosis, graphene oxide, plant
competition models or an open letter to Xi Jinping; those were replaced by
PubMed papers verified to mention the rsID.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts" / "setup"))

import audit_pmids as ap  # noqa: E402


@pytest.fixture(scope="module")
def audit_state():
    claims = ap.collect_claims(ap.DEFAULT_PANEL, ap.DEFAULT_RECS)
    cache = ap.load_cache(ap.DEFAULT_CACHE)
    allow = ap.load_allowlist(ap.DEFAULT_ALLOWLIST)
    return claims, cache, allow


def test_every_cited_pmid_is_in_the_committed_cache(audit_state):
    claims, cache, _ = audit_state
    missing = sorted({c["pmid"] for c in claims} - set(cache))
    assert not missing, ("PMIDs cited but absent from tests/fixtures/pubmed_cache.json - run "
                         f"scripts/setup/audit_pmids.py --update-cache and commit the cache: {missing}")


def test_no_citation_contradicts_its_pubmed_record(audit_state):
    claims, cache, allow = audit_state
    report = ap.audit(claims, cache, allow)
    assert report["n_claims"] > 250
    bad = [f"PMID {f['pmid']} [{f['source']}] {f['key']}: {f['title'][:70]!r}" for f in report["findings"]]
    assert not bad, "citations whose PubMed record mentions neither rsID, gene nor trait:\n" + "\n".join(bad)


def test_panel_rows_are_backed_by_gene_or_rsid_not_only_trait_words(audit_state):
    """Panel rows cite the association itself: at least 80% must match on
    rsID/gene, not merely on a trait word (weak). The weak remainder are
    consortium GWAS papers (GLGC lipids, GIANT BMI, MAGIC glycaemia, Dubois
    celiac) whose abstracts list no individual SNP."""
    claims, cache, allow = audit_state
    report = ap.audit(claims, cache, allow)
    panel = [r for r in report["results"] if r["source"] == "panel"]
    related = sum(1 for r in panel if r["status"] == "related")
    assert related / len(panel) >= 0.8, f"{related}/{len(panel)} panel citations name the rsID or gene"


def test_allowlist_entries_are_still_needed(audit_state):
    """An allowlisted citation that now passes on its own should be removed."""
    claims, cache, _ = audit_state
    report = ap.audit(claims, cache, {})
    flagged = {(f["pmid"], f["key"]) for f in report["findings"]}
    for key in ap.load_allowlist(ap.DEFAULT_ALLOWLIST):
        assert key in flagged, f"allowlist entry {key} is no longer flagged - delete it"


# ── unit tests of the matching rules ────────────────────────────────────────

def _claim(source="panel", tokens=("rs1", "FOXO3"), fallback=("longevity",)):
    return {"source": source, "key": "k", "pmid": "1", "tokens": list(tokens), "fallback_tokens": list(fallback)}


def test_gene_symbol_matches_with_suffix_and_alias():
    rec = {"title": "FOXO3A genotype is strongly associated with human longevity.", "abstract": ""}
    assert ap.evaluate(_claim(), rec)["status"] == "related"
    rec2 = {"title": "Apolipoprotein E polymorphism and atherosclerosis.", "abstract": ""}
    assert ap.evaluate(_claim(tokens=("rs429358", "APOE")), rec2)["status"] == "related"


def test_trait_words_only_give_weak_and_nothing_gives_unrelated():
    rec = {"title": "Genetics of human longevity: a review", "abstract": "no gene named"}
    assert ap.evaluate(_claim(), rec)["status"] == "weak"
    rec2 = {"title": "Rhinocerebral mucormycosis.", "abstract": "A fungal infection."}
    assert ap.evaluate(_claim(), rec2)["status"] == "unrelated"
    assert ap.evaluate(_claim(), None)["status"] == "missing"


def test_curation_history_pmids_are_not_treated_as_citations(tmp_path):
    import json
    panel = tmp_path / "panel.csv"
    panel.write_text("rsid,gene,trait_category,pmid\nrs1,GENE1,Trait one,111111\n")
    recs = tmp_path / "recs.json"
    recs.write_text(json.dumps({"Trait one": {
        "recommendation_en": "Your genotype (GENE1) ...", "recommendation_es": "",
        "reference": "Smith 2020, PMID:222222 (already in panel; corrected 2026-08-05 from a stale/wrong PMID 333333); "
                     "Doe 2021, PMID:444444 (corrected from a fabricated PMID:555555)"}}))
    claims = ap.collect_claims(panel, recs)
    rec_pmids = sorted(c["pmid"] for c in claims if c["source"] == "recommendation")
    assert rec_pmids == ["222222", "444444"]
    assert claims[0]["tokens"] == ["rs1", "GENE1"]
    assert "GENE1" in [c for c in claims if c["source"] == "recommendation"][0]["tokens"]
