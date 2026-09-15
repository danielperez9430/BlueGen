"""
Genotyping-array raw data → VCF (RELEASE_PLAN 3.0.3).

Accepts the consumer raw-data exports most people have — 23andMe,
AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA (23andMe-like) — and
writes a GRCh37 VCF the pipeline can run on. Unlike a WGS VCF, an array
file records EVERY assayed site, including homozygous-reference calls, so
those are emitted as 0/0 records: downstream stages then know which sites
were genotyped and the joint scoring (bluegen.joint, site_policy="array")
restricts itself to them instead of assuming absent == hom-ref.

Per row: chromosome codes are normalised (AncestryDNA 23/24/25/26 → X/Y/X/MT),
alleles must be A/C/G/T (no-calls '--'/'00' and indel codes I/D are counted
and skipped), the hg19 REF base is read from the FASTA and the genotype is
expressed against it; a row whose alleles contain no REF and no consistent
ALT is dropped as ref_mismatch. Haploid calls (one letter: X in males, Y,
MT) are written as homozygous. Files on GRCh38 are lifted row by row with
the hg38→hg19 chain (SNPs only) before the REF lookup.

Build detection: the file header ("build 37", "GRCh38", …) when present,
otherwise the positions of the curated-panel rsIDs found in the file are
compared with the panel's GRCh37 positions (and with their hg38 lift).
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from .genome_build import FastaIndex, GRCH37_LENGTHS, chain_chrom, normalize_chrom

FORMATS = ("23andme", "ancestrydna", "myheritage", "ftdna", "livingdna", "generic")
VALID = {"A", "C", "G", "T"}
_ANCESTRY_CHROM = {"23": "X", "24": "Y", "25": "X", "26": "MT"}
_BUILD37 = re.compile(r"build\s*37|GRCh37|hg19|b37", re.I)
_BUILD38 = re.compile(r"build\s*38|GRCh38|hg38|b38", re.I)


# ── reading ──────────────────────────────────────────────────────────────────

def open_text(path):
    """Plain, .gz or .zip (first member) → text handle."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", errors="replace")
    if p.suffix == ".zip":
        zf = zipfile.ZipFile(p)
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if not names:
            raise ValueError(f"{p}: empty zip")
        return io.TextIOWrapper(zf.open(names[0]), errors="replace")
    return open(p, errors="replace")


def sniff_format(path) -> tuple:
    """Returns (format, build_hint) from the first lines. build_hint ∈
    {'GRCh37','GRCh38',None}."""
    head = []
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            head.append(line.rstrip("\n"))
            if i > 60:
                break
    text = "\n".join(head)
    hint = "GRCh37" if _BUILD37.search(text) else ("GRCh38" if _BUILD38.search(text) else None)
    low = text.lower()
    if "23andme" in low:
        fmt = "23andme"
    elif "ancestrydna" in low or ("allele1" in low and "allele2" in low):
        fmt = "ancestrydna"
    elif "myheritage" in low:
        fmt = "myheritage"
    elif "familytreedna" in low or "ftdna" in low:
        fmt = "ftdna"
    elif "living dna" in low or "livingdna" in low:
        fmt = "livingdna"
    else:
        fmt = "generic"
    return fmt, hint


def _norm_chrom(raw: str) -> str:
    c = raw.strip().strip('"')
    c = _ANCESTRY_CHROM.get(c, c)
    return normalize_chrom(c)


def iter_genotypes(path, fmt=None):
    """Yields (rsid, chrom, pos, allele1, allele2) with chrom normalised
    ('1'..'22','X','Y','MT') and alleles upper-case strings ('-' for no-call,
    'I'/'D' for indel codes). Header/comment lines are skipped."""
    fmt = fmt or sniff_format(path)[0]
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            s = line.strip()
            if not s:
                continue
            if "," in s and "\t" not in s:
                parts = next(csv.reader([s]))
            else:
                parts = s.split("\t")
            parts = [p.strip().strip('"') for p in parts]
            if not parts or not parts[0].lower().startswith("rs") and not parts[0].lower().startswith("i"):
                # header rows ("rsid", "RSID", ...) or unknown identifiers
                if parts[0].lower() in ("rsid", "# rsid", "snp", "marker"):
                    continue
                if not re.match(r"^(rs|i)\d+$", parts[0], re.I):
                    continue
            try:
                if fmt == "ancestrydna" and len(parts) >= 5:
                    rsid, chrom, pos, a1, a2 = parts[0], parts[1], int(parts[2]), parts[3], parts[4]
                elif len(parts) >= 4:
                    rsid, chrom, pos, geno = parts[0], parts[1], int(parts[2]), parts[3]
                    geno = geno.replace(" ", "")
                    if len(geno) == 0:
                        a1 = a2 = "-"
                    elif len(geno) == 1:
                        a1 = a2 = geno   # haploid call → homozygous
                    else:
                        a1, a2 = geno[0], geno[1]
                else:
                    continue
            except ValueError:
                continue
            a1, a2 = a1.upper(), a2.upper()
            if a1 in ("0", "N"):
                a1 = "-"
            if a2 in ("0", "N"):
                a2 = "-"
            yield rsid, _norm_chrom(chrom), pos, a1, a2


# ── build detection ──────────────────────────────────────────────────────────

def detect_array_build(path, panel_positions: dict, chain_19to38=None, fmt=None, hint=None,
                       max_rows: int = 2_000_000) -> tuple:
    """panel_positions: {rsid: (chrom, pos_grch37)}. Returns (build, evidence).
    The header hint wins when the marker check agrees or cannot decide;
    a strong marker disagreement overrides it (headers are copied around)."""
    if fmt is None or hint is None:
        fmt, hint = sniff_format(path)
    found = {}
    for i, (rsid, chrom, pos, _a1, _a2) in enumerate(iter_genotypes(path, fmt)):
        if rsid in panel_positions:
            found[rsid] = (chrom, pos)
            if len(found) >= len(panel_positions):
                break
        if i >= max_rows:
            break
    n = len(found)
    match37 = sum(1 for r, (c, p) in found.items()
                  if normalize_chrom(str(panel_positions[r][0])) == c and int(panel_positions[r][1]) == p)
    match38 = None
    if chain_19to38 is not None and n:
        from pyliftover import LiftOver
        lo = LiftOver(str(chain_19to38))
        match38 = 0
        for r, (c, p) in found.items():
            pc, pp = panel_positions[r]
            conv = lo.convert_coordinate(chain_chrom(normalize_chrom(str(pc))), int(pp) - 1)
            if conv and normalize_chrom(conv[0][0]) == c and conv[0][1] + 1 == p:
                match38 += 1
    evidence = {"format": fmt, "header_hint": hint, "panel_markers_found": n,
                "match_grch37": match37, "match_grch38": match38}
    if n >= 10:
        frac37 = match37 / n
        frac38 = (match38 / n) if match38 is not None else None
        if frac37 >= 0.9 and (frac38 is None or frac37 > frac38):
            evidence["method"] = "marker_positions"
            return "GRCh37", evidence
        if frac38 is not None and frac38 >= 0.9 and frac38 > frac37:
            evidence["method"] = "marker_positions"
            return "GRCh38", evidence
    if hint:
        evidence["method"] = "header_hint"
        return hint, evidence
    evidence["method"] = None
    return "unknown", evidence


# ── conversion ───────────────────────────────────────────────────────────────

def array_to_vcf(path, out_vcf, fasta_path, sample_id: str, fmt=None, lifter=None,
                 bcftools: str = "bcftools", progress=None) -> dict:
    """Write a GRCh37 VCF (bgzipped, sorted, indexed) from an array file.
    `lifter` is a pyliftover.LiftOver(hg38ToHg19) for GRCh38 files. Every
    assayed site is written, homozygous-reference ones as 0/0."""
    fmt = fmt or sniff_format(path)[0]
    fasta = FastaIndex(fasta_path)
    out_vcf = Path(out_vcf)
    stats = {k: 0 for k in ("n_rows", "n_written", "n_hom_ref", "n_het", "n_hom_alt", "n_multiallelic",
                            "n_nocall", "n_indel_code", "n_ref_mismatch", "n_unknown_chrom",
                            "n_unmapped", "n_strand_flipped", "n_duplicate_site")}
    stats["format"] = fmt
    stats["lifted_from_grch38"] = lifter is not None

    tmp_dir = Path(tempfile.mkdtemp(prefix="bluegen_array_"))
    tmp_vcf = tmp_dir / "array_unsorted.vcf"
    seen = set()
    comp = str.maketrans("ACGT", "TGCA")
    with open(tmp_vcf, "w") as dst:
        dst.write("##fileformat=VCFv4.2\n")
        dst.write("##reference=GRCh37 (hg19)\n")
        dst.write(f"##source=BlueGen array_to_vcf format={fmt} lifted_from_grch38={lifter is not None}\n")
        dst.write("##INFO=<ID=RSID,Number=1,Type=String,Description=\"Array marker id\">\n")
        dst.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for c, length in GRCH37_LENGTHS.items():
            dst.write(f"##contig=<ID={c},length={length}>\n")
        dst.write("##contig=<ID=MT,length=16569>\n")
        dst.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")

        for rsid, chrom, pos, a1, a2 in iter_genotypes(path, fmt):
            stats["n_rows"] += 1
            if progress and stats["n_rows"] % 200000 == 0:
                progress(stats["n_rows"])
            if chrom not in GRCH37_LENGTHS and chrom != "MT":
                stats["n_unknown_chrom"] += 1
                continue
            if a1 == "-" or a2 == "-":
                stats["n_nocall"] += 1
                continue
            if a1 in ("I", "D") or a2 in ("I", "D"):
                stats["n_indel_code"] += 1
                continue
            if a1 not in VALID or a2 not in VALID:
                stats["n_nocall"] += 1
                continue
            if lifter is not None:
                conv = lifter.convert_coordinate(chain_chrom(chrom), pos - 1)
                if len(conv) != 1 or normalize_chrom(conv[0][0]) != chrom:
                    stats["n_unmapped"] += 1
                    continue
                pos = conv[0][1] + 1
                if conv[0][2] == "-":
                    a1, a2 = a1.translate(comp), a2.translate(comp)
                    stats["n_strand_flipped"] += 1
            key = (chrom, pos)
            if key in seen:
                stats["n_duplicate_site"] += 1
                continue
            fname = chain_chrom(chrom)
            if not fasta.has(fname):
                stats["n_unknown_chrom"] += 1
                continue
            ref = fasta.fetch(fname, pos, 1)
            if ref not in VALID:
                stats["n_ref_mismatch"] += 1
                continue
            alts = []
            gt_idx = []
            for a in (a1, a2):
                if a == ref:
                    gt_idx.append(0)
                else:
                    if a not in alts:
                        alts.append(a)
                    gt_idx.append(alts.index(a) + 1)
            if len(alts) == 2:
                stats["n_multiallelic"] += 1
            if gt_idx == [0, 0]:
                stats["n_hom_ref"] += 1
            elif 0 in gt_idx:
                stats["n_het"] += 1
            else:
                stats["n_hom_alt"] += 1 if len(alts) == 1 else 0
            seen.add(key)
            alt_field = ",".join(alts) if alts else "."
            dst.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt_field}\t.\tPASS\tRSID={rsid}\tGT\t"
                      f"{gt_idx[0]}/{gt_idx[1]}\n")
            stats["n_written"] += 1
    fasta.close()

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    sort_tmp = tmp_dir / "sort"
    sort_tmp.mkdir()
    r = subprocess.run([bcftools, "sort", "-T", str(sort_tmp), "-Oz", "-o", str(out_vcf), str(tmp_vcf)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise RuntimeError(f"bcftools sort failed: {r.stderr[-400:]}")
    r = subprocess.run([bcftools, "index", "-f", "-t", str(out_vcf)], capture_output=True, text=True)
    if r.returncode != 0:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise RuntimeError(f"bcftools index failed: {r.stderr[-400:]}")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    stats["output"] = str(out_vcf)
    return stats


def panel_positions_from_csv(snp_db_path) -> dict:
    """{rsid: (chrom, pos)} from the curated panel CSV (GRCh37)."""
    out = {}
    with open(snp_db_path) as fh:
        for row in csv.DictReader(fh):
            rsid = (row.get("rsid") or "").strip()
            chrom = (row.get("chrom") or "").strip()
            pos = (row.get("pos") or "").strip()
            try:
                if rsid and chrom and pos and float(pos) > 0:
                    out[rsid] = (chrom, int(float(pos)))
            except ValueError:
                continue
    return out


__all__ = ["FORMATS", "sniff_format", "iter_genotypes", "detect_array_build", "array_to_vcf",
           "panel_positions_from_csv"]
