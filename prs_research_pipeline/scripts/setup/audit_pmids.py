#!/usr/bin/env python3
"""
Audit the literature citations of the curated panel and the trait
recommendations (RELEASE_PLAN 3.0.6).

History: rounds 9-11 of IMPROVEMENT_PLAN 1.2.C / 1.4 found that wrong PMIDs
were SYSTEMIC in the panel - rows citing an unrelated crystallography or
education-methodology paper - and every one was caught by hand. This script
makes that check mechanical: for every PMID cited by a panel row
(data/snp_database_annotated.csv) or a recommendation
(data/trait_recommendations.json), fetch title + abstract from PubMed
(E-utilities efetch) and verify the text mentions what the citation is for:
the rsID, one of the gene symbols, or - for recommendations without a gene
in the text - the trait's key words.

A JSON cache (default tests/fixtures/pubmed_cache.json, committed) makes the
check reproducible offline; tests/test_pmid_audit.py runs the same logic
against it in CI without any network call. Justified exceptions (a GWAS of
the trait whose abstract names no gene) live in
data/pmid_audit_allowlist.json with a reason each.

Usage:
    python audit_pmids.py                      # offline against the cache, exit 1 on findings
    python audit_pmids.py --update-cache       # fetch missing PMIDs from PubMed, then audit
    python audit_pmids.py --json report.json   # + machine-readable report
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[2]
REPO = PIPELINE.parent
DEFAULT_PANEL = PIPELINE / "data" / "snp_database_annotated.csv"
DEFAULT_RECS = PIPELINE / "data" / "trait_recommendations.json"
DEFAULT_CACHE = REPO / "tests" / "fixtures" / "pubmed_cache.json"
DEFAULT_ALLOWLIST = PIPELINE / "data" / "pmid_audit_allowlist.json"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 100
STOPWORDS = {"levels", "status", "risk", "response", "metabolism", "sensitivity", "function", "type",
             "extended", "early", "male", "pattern", "intolerance", "predisposition", "regulation",
             "performance", "reaction", "adjusted", "body", "index", "mass", "cognition", "neuro"}
# Older papers write the protein/enzyme name instead of the HGNC symbol.
GENE_ALIASES = {
    "APOE": ["apolipoprotein E", "apo E", "apoE"], "APOA1": ["apolipoprotein A-I", "apo A-I", "apoA-I"],
    "APOA5": ["apolipoprotein A-V", "apoA-V"], "APOB": ["apolipoprotein B", "apo B", "apoB"],
    "LPL": ["lipoprotein lipase"], "LDLR": ["LDL receptor", "low density lipoprotein receptor",
                                             "low-density lipoprotein receptor"],
    "PCSK9": ["proprotein convertase"], "MTHFR": ["methylenetetrahydrofolate reductase"],
    "MTR": ["methionine synthase"], "MTRR": ["methionine synthase reductase"],
    "CYP1A2": ["cytochrome P450 1A2", "CYP1A2"], "CYP2D6": ["cytochrome P450 2D6"],
    "CYP2C19": ["cytochrome P450 2C19"], "CYP2C9": ["cytochrome P450 2C9"],
    "VKORC1": ["vitamin K epoxide reductase"], "ADH1B": ["alcohol dehydrogenase"],
    "ALDH2": ["aldehyde dehydrogenase"], "LCT": ["lactase"], "MCM6": ["lactase"],
    "FTO": ["fat mass and obesity"], "MC4R": ["melanocortin-4 receptor", "melanocortin 4 receptor"],
    "MC1R": ["melanocortin-1 receptor", "melanocortin 1 receptor", "redhead", "red hair"],
    "TCF7L2": ["transcription factor 7-like 2"], "PPARG": ["PPAR-gamma", "PPARgamma", "peroxisome proliferator"],
    "ACE": ["angiotensin-converting enzyme", "angiotensin converting enzyme"],
    "ACTN3": ["alpha-actinin-3", "α-actinin-3", "actinin"], "COMT": ["catechol-O-methyltransferase"],
    "BDNF": ["brain-derived neurotrophic factor"], "OPRM1": ["mu-opioid receptor", "μ-opioid receptor"],
    "ADORA2A": ["adenosine A2A receptor", "A2a receptor"], "HLA-DQA1": ["HLA-DQ"], "HLA-DQB1": ["HLA-DQ"],
    "HFE": ["hemochromatosis", "haemochromatosis"], "TMPRSS6": ["matriptase-2"],
    "SLC23A1": ["sodium-dependent vitamin C transporter", "SVCT1"], "GC": ["vitamin D binding protein", "vitamin D-binding protein"],
    "CYP2R1": ["25-hydroxylase"], "VDR": ["vitamin D receptor"], "FADS1": ["fatty acid desaturase"],
    "FADS2": ["fatty acid desaturase"], "TNF": ["tumor necrosis factor", "tumour necrosis factor", "TNF-alpha", "TNFα"],
    "IL6": ["interleukin-6", "interleukin 6", "IL-6"], "IL1B": ["interleukin-1", "IL-1"], "CRP": ["C-reactive protein"],
    "NOS3": ["endothelial nitric oxide synthase", "eNOS"], "AGT": ["angiotensinogen"],
    "SLC30A8": ["zinc transporter"], "GCKR": ["glucokinase regulatory"], "G6PC2": ["glucose-6-phosphatase"],
    "ABCA1": ["ATP-binding cassette"], "CETP": ["cholesteryl ester transfer protein"],
    "AR": ["androgen receptor"], "EDA2R": ["ectodysplasin"], "TERC": ["telomerase RNA"], "TERT": ["telomerase reverse transcriptase"],
    "OBFC1": ["STN1"], "FOXO3": ["FOXO3A", "forkhead"], "CLOCK": ["circadian"], "PER2": ["period 2", "Period2"],
    "HNMT": ["histamine N-methyltransferase"], "AOC1": ["diamine oxidase", "DAO"], "ABO": ["blood group"],
    "FUT2": ["secretor"], "SLC45A2": ["MATP"], "OCA2": ["HERC2", "P protein"], "HERC2": ["OCA2"],
    "IRF4": ["interferon regulatory factor 4"], "TYR": ["tyrosinase"], "SLC24A5": ["NCKX5"],
    "ABCC11": ["earwax", "ear wax", "cerumen"], "GSTP1": ["glutathione S-transferase"], "GSTM1": ["glutathione S-transferase"],
    "GSTT1": ["glutathione S-transferase"], "NAT2": ["N-acetyltransferase"], "CYP1A1": ["cytochrome P450 1A1"],
    "SOD2": ["manganese superoxide dismutase", "MnSOD"], "TAS2R38": ["PTC", "phenylthiocarbamide", "bitter taste"],
    "OXTR": ["oxytocin receptor"], "DRD2": ["dopamine D2 receptor", "D2 receptor"], "DRD4": ["dopamine D4 receptor"],
    "ANKK1": ["Taq1A", "TaqIA"], "SLC6A4": ["serotonin transporter", "5-HTTLPR"], "MAOA": ["monoamine oxidase"],
    "ADRB2": ["beta2-adrenergic receptor", "β2-adrenergic"], "ADRB3": ["beta3-adrenergic receptor", "β3-adrenergic"],
    "UCP1": ["uncoupling protein"], "LEPR": ["leptin receptor"], "LEP": ["leptin"], "ADIPOQ": ["adiponectin"],
    "KIBRA": ["WWC1"], "WWC1": ["KIBRA"], "APOA2": ["apolipoprotein A-II"], "LIPC": ["hepatic lipase"],
    "SLC2A2": ["GLUT2"], "TAS1R2": ["sweet taste"], "CD36": ["fatty acid translocase"], "PEMT": ["phosphatidylethanolamine"],
    "CHDH": ["choline dehydrogenase"], "BCMO1": ["beta-carotene", "β-carotene", "BCO1"], "BCO1": ["beta-carotene", "BCMO1"],
    "TCN2": ["transcobalamin"], "FUT6": ["fucosyltransferase"], "CUBN": ["cubilin"], "TF": ["transferrin"],
    "SLC17A1": ["urate"], "SLC2A9": ["urate", "GLUT9"], "ABCG2": ["urate", "BCRP"], "COL1A1": ["collagen"],
    "MMP1": ["matrix metalloproteinase"], "MMP3": ["matrix metalloproteinase"], "CLEC4E": ["Mincle"],
    "SLC30A2": ["zinc transporter"], "SLC39A8": ["zinc transporter", "ZIP8"], "CBS": ["cystathionine"],
    "SLC19A1": ["reduced folate carrier"], "DHFR": ["dihydrofolate reductase"], "AHR": ["aryl hydrocarbon receptor"],
    "RHD": ["Rh", "RhD"], "SHBG": ["sex hormone-binding globulin"], "CYP19A1": ["aromatase"],
    "ESR1": ["estrogen receptor", "oestrogen receptor"], "HSD17B1": ["17beta-hydroxysteroid"],
    "SLC4A5": ["salt sensitivity"], "ADD1": ["adducin"], "CYP11B2": ["aldosterone synthase"],
}
# Trait names → extra key words a paper about the trait would use.
TRAIT_SYNONYMS = {
    "urate levels": ["gout", "uric acid", "urate"], "tendinopathy risk": ["tendon", "achilles", "tendinopathy", "collagen"],
    "alzheimer risk": ["dementia", "cognitive", "alzheimer"], "apolipoprotein b": ["LDL", "hypercholesterolemia", "apoB", "cholesterol"],
    "baldness (male pattern)": ["alopecia", "hair loss", "androgenetic"], "hair color (red)": ["redhead", "red hair", "MC1R"],
    "blood pressure": ["hypertension", "blood pressure"], "telomere length": ["telomer"],
    "alcohol flush reaction": ["alcohol", "acetaldehyde", "flush"], "vitamin d metabolism": ["vitamin D", "25-hydroxyvitamin", "cholecalciferol"],
    "cognitive function": ["cognitive", "memory", "hippocamp"], "inflammation (crp levels)": ["C-reactive", "CRP", "inflammat"],
    "caffeine metabolism": ["caffeine", "coffee"], "caffeine sensitivity (anxiety)": ["caffeine", "anxiety"],
    "lactose intolerance": ["lactose", "lactase", "hypolactasia"], "folate & methylation": ["folate", "homocysteine", "methyl"],
    "glucose metabolism": ["glucose", "diabetes", "insulin"], "lipid metabolism": ["lipid", "cholesterol", "triglycerid"],
    "hdl/ldl extended": ["lipid", "cholesterol", "HDL", "LDL"], "omega-3 metabolism": ["omega-3", "fatty acid", "PUFA"],
    "detoxification": ["xenobiotic", "detoxif", "glutathione", "carcinogen"], "iron levels": ["iron", "ferritin", "hepcidin"],
    "histamine intolerance": ["histamine", "diamine oxidase"], "muscle performance": ["muscle", "athlet", "endurance", "sprint"],
    "obesity predisposition": ["obesity", "BMI", "body mass"], "vitamin c": ["ascorb", "vitamin C"],
}
_PMID_RE = re.compile(r"PMID:?\s*(\d{6,9})")
_HISTORY_RE = re.compile(r"(?:wrong|stale|fabricated|previous|corrected(?:\s+\S+){0,6}?\s+from)[^;)]{0,80}?PMID:?\s*\d{6,9}", re.I)
_GENE_PAREN_RE = re.compile(r"\(([A-Z][A-Z0-9]{1,9}(?:\s*/\s*[A-Z][A-Z0-9]{1,9})*)\)")
_RSID_RE = re.compile(r"\brs\d+\b")


# ── claims ───────────────────────────────────────────────────────────────────

def _split_genes(gene_field: str) -> list:
    return [g.strip() for g in re.split(r"[/ ,;]+", gene_field or "") if len(g.strip()) >= 3]


def _trait_keywords(trait: str) -> list:
    words = re.findall(r"[A-Za-z][A-Za-z-]{3,}", trait or "")
    out = [w for w in words if w.lower() not in STOPWORDS]
    return out + TRAIT_SYNONYMS.get((trait or "").strip().lower(), [])


def collect_claims(panel_csv=DEFAULT_PANEL, recs_json=DEFAULT_RECS) -> list:
    """Each claim: {source, key, pmid, tokens, fallback_tokens}. `tokens` are
    the strong evidence (rsID, gene symbols); `fallback_tokens` (trait words)
    only count when the claim carries no gene at all."""
    claims = []
    trait_tokens = {}   # trait_category (lower) → panel rsIDs + genes, for the recommendations
    with open(panel_csv) as fh:
        for row in csv.DictReader(fh):
            tkey = (row.get("trait_category") or "").strip().lower()
            toks = [row["rsid"]] + _split_genes(row.get("gene", ""))
            trait_tokens.setdefault(tkey, [])
            trait_tokens[tkey] += [t for t in toks if t not in trait_tokens[tkey]]
            pmid = (row.get("pmid") or "").strip()
            if not pmid or not pmid.isdigit():
                continue
            claims.append({
                "source": "panel", "key": f"{row['rsid']} ({row.get('gene', '')}, {row.get('trait_category', '')})",
                "pmid": pmid,
                "tokens": toks,
                "fallback_tokens": _trait_keywords(row.get("trait_category", "")),
            })
    recs = json.loads(Path(recs_json).read_text())
    for trait, rec in recs.items():
        if trait.startswith("_") or not isinstance(rec, dict):
            continue
        text = f"{rec.get('recommendation_en', '')} {rec.get('recommendation_es', '')}"
        genes = []
        for m in _GENE_PAREN_RE.finditer(text):
            genes += _split_genes(m.group(1))
        rsids = sorted(set(_RSID_RE.findall(text)))
        tokens = rsids + genes + trait_tokens.get(trait.strip().lower(), [])
        tokens = list(dict.fromkeys(tokens))
        # The reference prose keeps curation history ("corrected from a
        # fabricated PMID:12771985"): those PMIDs are not citations.
        ref_text = _HISTORY_RE.sub(" ", rec.get("reference", ""))
        for pmid in sorted(set(_PMID_RE.findall(ref_text))):
            claims.append({
                "source": "recommendation", "key": trait, "pmid": pmid,
                "tokens": tokens,
                "fallback_tokens": _trait_keywords(trait),
            })
    return claims


# ── PubMed ───────────────────────────────────────────────────────────────────

def fetch_pubmed(pmids: list, sleep: float = 0.4) -> dict:
    out = {}
    pmids = sorted(set(pmids))
    for i in range(0, len(pmids), BATCH):
        chunk = pmids[i:i + BATCH]
        url = f"{EFETCH}?db=pubmed&id={','.join(chunk)}&rettype=abstract&retmode=xml"
        req = urllib.request.Request(url, headers={"User-Agent": "BlueGen-pmid-audit/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            root = ET.fromstring(r.read())
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//MedlineCitation/PMID") or art.findtext(".//PMID")
            title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
            abstract = " ".join("".join(t.itertext()) for t in art.findall(".//AbstractText"))
            year = art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate") or ""
            out[pmid] = {"title": title.strip(), "abstract": abstract.strip()[:4000], "year": year}
        if i + BATCH < len(pmids):
            time.sleep(sleep)
    return out


ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def find_paper_for_rsid(rsid: str, genes: list, cache: dict, trait_words: list = (), retmax: int = 40) -> tuple:
    """PubMed search for the rsID; returns (pmid, record, why) for the best
    candidate whose title/abstract actually contains the rsID, ranked:
      tier 3  rsID + a gene of the row + a trait key word (the association itself)
      tier 2  rsID + gene
      tier 1  rsID only
    Within a tier the oldest paper wins (primary reports predate reviews and
    meta-analyses). (None, None, why) when nothing mentions the rsID."""
    from urllib.parse import quote
    url = f"{ESEARCH}?db=pubmed&term={quote(rsid + '[All Fields]')}&retmode=json&retmax={retmax}&sort=relevance"
    req = urllib.request.Request(url, headers={"User-Agent": "BlueGen-pmid-audit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        ids = json.loads(r.read()).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None, None, "no PubMed record mentions this rsID"
    missing = [i for i in ids if i not in cache]
    if missing:
        cache.update(fetch_pubmed(missing))
        time.sleep(0.4)

    def year_of(rec):
        m = re.search(r"\d{4}", rec.get("year", "") or "")
        return int(m.group()) if m else 9999

    candidates = []
    for pmid in ids:
        rec = cache.get(pmid)
        if not rec:
            continue
        text = f"{rec['title']} {rec['abstract']}"
        if not _mentions(text, rsid):
            continue
        gene_hit = any(_mentions(text, g) or any(_mentions(text, a) for a in GENE_ALIASES.get(g.upper(), []))
                       for g in genes)
        trait_hit = any(_mentions(text, w) for w in trait_words)
        tier = 3 if (gene_hit and trait_hit) else (2 if gene_hit else 1)
        candidates.append((tier, -year_of(rec), pmid, rec))
    if not candidates:
        return None, None, "PubMed hits for the rsID do not mention it in title/abstract"
    tier, _neg_year, pmid, rec = max(candidates, key=lambda c: (c[0], c[1]))
    why = {3: "mentions the rsID, the gene and the trait", 2: "mentions the rsID and the gene",
           1: "mentions the rsID"}[tier]
    return pmid, rec, why


def load_cache(path) -> dict:
    p = Path(path)
    return json.loads(p.read_text()) if p.exists() else {}


def save_cache(path, cache: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(sorted(cache.items())), indent=1, ensure_ascii=False) + "\n")


# ── evaluation ───────────────────────────────────────────────────────────────

def _mentions(text: str, token: str) -> bool:
    if token.lower().startswith("rs") and token[2:].isdigit():
        return re.search(rf"\b{re.escape(token)}\b", text, re.I) is not None
    # gene symbols may carry a suffix in the literature (FOXO3A, CYP2D6*4,
    # ADORA2A); require a clean start and at most one trailing letter/star
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(token)}[A-Za-z*]?(?![A-Za-z0-9])", text, re.I) is not None


def evaluate(claim: dict, record: dict | None) -> dict:
    """status: related | weak | unrelated | missing.

    related   the record names the rsID or a gene of the claim;
    weak      only trait key words matched - accepted for recommendations,
              whose citations often back the intervention rather than the
              gene (a Mediterranean-diet/telomere trial names no TERC), and
              for panel rows without any gene/rsID token;
    unrelated nothing matched - the finding this audit exists for."""
    if not record:
        return {**claim, "status": "missing", "title": ""}
    text = f"{record.get('title', '')} {record.get('abstract', '')}"
    hits = [t for t in claim["tokens"] if _mentions(text, t)]
    if hits:
        return {**claim, "status": "related", "matched": hits, "title": record.get("title", "")}
    alias_hits = [a for t in claim["tokens"] for a in GENE_ALIASES.get(t.upper(), []) if _mentions(text, a)]
    if alias_hits:
        return {**claim, "status": "related", "matched": alias_hits, "title": record.get("title", "")}
    fb = [t for t in claim["fallback_tokens"] if _mentions(text, t)]
    if fb:
        return {**claim, "status": "weak", "matched": fb, "title": record.get("title", "")}
    return {**claim, "status": "unrelated", "matched": fb, "title": record.get("title", "")}


def load_allowlist(path) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    entries = json.loads(p.read_text())
    return {(e["pmid"], e["key"]): e.get("reason", "") for e in entries}


def audit(claims: list, cache: dict, allowlist: dict) -> dict:
    results = [evaluate(c, cache.get(c["pmid"])) for c in claims]
    findings = []
    for r in results:
        if r["status"] in ("unrelated", "missing"):
            reason = allowlist.get((r["pmid"], r["key"]))
            if reason:
                r["status"] = "allowlisted"
                r["reason"] = reason
            else:
                findings.append(r)
    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    return {"n_claims": len(results), "n_pmids": len({r["pmid"] for r in results}), "counts": counts,
            "findings": findings, "results": results}


def fix_panel_pmids(report: dict, cache: dict, panel_csv, force_rsids=()) -> list:
    """Replace or blank the PMID of every panel row flagged unrelated/missing
    (plus any row of `force_rsids`, e.g. a 'weak' match judged wrong by hand),
    in every data/snp_database*.csv that carries the same (rsid, trait) row."""
    targets = {}
    force = set(force_rsids)
    rows = [f for f in report["findings"] if f["source"] == "panel"]
    rows += [r for r in report["results"] if r["source"] == "panel" and r["tokens"][0] in force
             and r["status"] not in ("unrelated", "missing")]
    for f in rows:
        rsid = f["tokens"][0]
        genes = f["tokens"][1:]
        key = (rsid, f["pmid"])
        if key in targets:
            continue
        new_pmid, rec, why = find_paper_for_rsid(rsid, genes, cache, trait_words=f.get("fallback_tokens", []))
        if rsid in force and new_pmid == f["pmid"]:
            continue   # the current paper is already the best candidate
        targets[key] = {"rsid": rsid, "gene": "/".join(genes), "old_pmid": f["pmid"], "old_title": f["title"],
                        "new_pmid": new_pmid or "", "new_title": rec["title"] if rec else "", "why": why}
        time.sleep(0.35)
    if not targets:
        return []
    data_dir = Path(panel_csv).parent
    for csv_path in sorted(data_dir.glob("snp_database*.csv")):
        # snp_database.csv opens with '#' comment lines before the header:
        # keep that preamble byte-for-byte, rewrite only the table.
        raw = csv_path.read_text(newline="").splitlines(keepends=True)
        hdr_idx = next((i for i, ln in enumerate(raw) if ln.startswith("rsid,")), None)
        if hdr_idx is None:
            continue
        preamble = "".join(raw[:hdr_idx])
        import io
        reader = csv.DictReader(io.StringIO("".join(raw[hdr_idx:])))
        fieldnames = reader.fieldnames or []
        rows = list(reader)
        if "pmid" not in fieldnames:
            continue
        changed = 0
        for row in rows:
            key = (row.get("rsid", ""), (row.get("pmid") or "").strip())
            if key in targets:
                row["pmid"] = targets[key]["new_pmid"]
                changed += 1
        if changed:
            bak = csv_path.with_suffix(csv_path.suffix + ".bak")
            bak.write_bytes(csv_path.read_bytes())
            with open(csv_path, "w", newline="") as fh:
                fh.write(preamble)
                w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            print(f"  ✎ {csv_path.name}: {changed} row(s) updated (backup: {bak.name})")
    return list(targets.values())


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    ap.add_argument("--recommendations", type=Path, default=DEFAULT_RECS)
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--allowlist", type=Path, default=DEFAULT_ALLOWLIST)
    ap.add_argument("--update-cache", action="store_true", help="fetch PMIDs missing from the cache")
    ap.add_argument("--json", type=Path, help="write the machine-readable report here")
    ap.add_argument("--fix", action="store_true",
                    help="panel rows with an unrelated/missing PMID: replace it with a PubMed paper that mentions "
                         "the rsID (relevance-ranked, gene preferred) or blank it; rewrites every data/snp_database*.csv "
                         "carrying the row (keeps .bak) and writes <json>.fixes.json")
    ap.add_argument("--fix-rsids", default="",
                    help="comma-separated rsIDs to re-source even if not flagged (weak matches judged wrong by hand)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    claims = collect_claims(args.panel, args.recommendations)
    cache = load_cache(args.cache)
    missing = sorted({c["pmid"] for c in claims} - set(cache))
    if missing and args.update_cache:
        print(f"Fetching {len(missing)} PMIDs from PubMed…")
        cache.update(fetch_pubmed(missing))
        save_cache(args.cache, cache)
        still = sorted({c["pmid"] for c in claims} - set(cache))
        if still:
            print(f"  ⚠️  {len(still)} PMIDs returned nothing (invalid?): {', '.join(still[:10])}")
    elif missing:
        print(f"⚠️  {len(missing)} PMIDs not in cache ({args.cache}); run with --update-cache")

    report = audit(claims, cache, load_allowlist(args.allowlist))

    if args.fix:
        fixes = fix_panel_pmids(report, cache, args.panel,
                                force_rsids=[r.strip() for r in args.fix_rsids.split(",") if r.strip()])
        save_cache(args.cache, cache)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.with_suffix(".fixes.json").write_text(json.dumps(fixes, indent=2, ensure_ascii=False))
        print(f"\n🔧 {len(fixes)} panel citation(s) changed "
              f"({sum(1 for f in fixes if f['new_pmid'])} replaced, {sum(1 for f in fixes if not f['new_pmid'])} blanked)")
        for f in fixes:
            print(f"  {f['rsid']:<12} {f['gene']:<9} {f['old_pmid']:<9} → {f['new_pmid'] or '(blank)':<9} {f['why']}"
                  + (f"\n               {f['new_title'][:100]}" if f.get("new_title") else ""))
        # re-audit after the rewrite
        claims = collect_claims(args.panel, args.recommendations)
        report = audit(claims, cache, load_allowlist(args.allowlist))

    if not args.quiet:
        print(f"\nPMID audit: {report['n_claims']} citations, {report['n_pmids']} unique PMIDs")
        for k, v in sorted(report["counts"].items()):
            print(f"  {k:<12}{v}")
        if report["findings"]:
            print(f"\n❌ {len(report['findings'])} citation(s) whose PubMed record mentions neither the rsID, "
                  "the gene nor the trait:")
            for f in report["findings"]:
                print(f"  PMID {f['pmid']:<9} [{f['source']}] {f['key'][:60]}")
                print(f"           expected: {', '.join(f['tokens'] or f['fallback_tokens'])}")
                print(f"           title:    {f['title'][:110] if f['title'] else '(no PubMed record)'}")
        else:
            print("\n✅ every citation is consistent with its PubMed record (or allowlisted with a reason)")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, ensure_ascii=False))
        print(f"📁 {args.json}")
    return 1 if report["findings"] else 0


if __name__ == "__main__":
    sys.exit(main())
