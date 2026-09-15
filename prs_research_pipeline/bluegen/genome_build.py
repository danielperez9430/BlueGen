"""
Genome-build detection and GRCh38 → GRCh37 liftover for input VCFs
(RELEASE_PLAN 3.0.2).

Every downstream stage of BlueGen is GRCh37/hg19: the curated panel
positions, the 1000 Genomes PLINK reference, ClinVar, PharmGKB, the PGS
Catalog harmonized files and the hg19 FASTA. A GRCh38 VCF (what Nebula,
Dante, Sequencing.com and most 2024+ pipelines deliver) used to enter
Stage A silently and score almost nothing. This module:

  1. `detect_build(vcf)` — decides GRCh37 / GRCh38 / unknown from the VCF
     header: contig lengths (strongest signal, present in every DeepVariant /
     GATK / bcftools header) and the ##reference line; an optional marker
     probe (`detect_build_from_markers`) resolves headerless files.
  2. `liftover_vcf(...)` — streams a GRCh38 VCF to GRCh37 through a UCSC
     chain file with pyliftover (pure Python), complementing alleles on
     minus-strand blocks, verifying REF against the hg19 FASTA, swapping
     REF/ALT (and flipping GT/AD/PL) where hg38 REF is the hg19 ALT, and
     dropping what cannot be mapped safely. Output is bgzipped, sorted and
     indexed with bcftools. Every drop is counted and reported.

Nothing here touches the pipeline's own outputs; prs.py decides what to
do with the result and records it in reproducibility/input_build.json.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path

# ── Reference contig lengths (autosomes + X/Y) ──────────────────────────────
GRCH37_LENGTHS = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}
GRCH38_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
    "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
    "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}
STANDARD_CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]

CHAIN_URLS = {
    "hg38ToHg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    "hg19ToHg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
}

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# ── chromosome naming ────────────────────────────────────────────────────────

def normalize_chrom(name: str) -> str:
    """'chr1' → '1', 'chrM'/'M' → 'MT', 'chrX' → 'X'. Unknown contigs unchanged."""
    n = name[3:] if name.lower().startswith("chr") else name
    if n in ("M", "MT"):
        return "MT"
    return n


def chain_chrom(norm: str) -> str:
    """UCSC chain/FASTA naming for a normalised chromosome ('1' → 'chr1', 'MT' → 'chrM')."""
    return "chrM" if norm == "MT" else f"chr{norm}"


def styled_chrom(norm: str, chr_prefix: bool) -> str:
    if not chr_prefix:
        return norm
    return "chrM" if norm == "MT" else f"chr{norm}"


# ── build detection ──────────────────────────────────────────────────────────

_REF37 = re.compile(r"hs37d5|GRCh37|hg19|\bb37\b|NCBI37|human_g1k_v37", re.I)
_REF38 = re.compile(r"GRCh38|hg38|\bb38\b|NCBI38", re.I)
_CONTIG = re.compile(r"^##contig=<(.*)>", re.M)


def detect_build_from_header(header: str) -> tuple:
    """Returns (build, evidence). build ∈ {'GRCh37', 'GRCh38', 'unknown'}.

    Contig lengths decide when at least 3 standard contigs carry a length:
    every matching length is a vote (a mixed header is 'unknown'). The
    ##reference / ##assembly text is used when no contig lengths exist."""
    votes = {"GRCh37": 0, "GRCh38": 0}
    n_contigs = 0
    for m in _CONTIG.finditer(header):
        fields = dict(kv.split("=", 1) for kv in m.group(1).split(",") if "=" in kv)
        cid = normalize_chrom(fields.get("ID", ""))
        if cid not in GRCH37_LENGTHS or "length" not in fields:
            continue
        try:
            length = int(fields["length"])
        except ValueError:
            continue
        n_contigs += 1
        if length == GRCH37_LENGTHS[cid]:
            votes["GRCh37"] += 1
        elif length == GRCH38_LENGTHS[cid]:
            votes["GRCh38"] += 1

    ref_lines = [ln for ln in header.splitlines() if ln.startswith(("##reference", "##assembly"))]
    ref_text = " ".join(ref_lines)
    ref_hint = "GRCh37" if _REF37.search(ref_text) else ("GRCh38" if _REF38.search(ref_text) else None)

    evidence = {"contigs_with_length": n_contigs, "votes": votes, "reference_line_hint": ref_hint,
                "method": None}
    if n_contigs >= 3 and votes["GRCh37"] != votes["GRCh38"]:
        winner = "GRCh37" if votes["GRCh37"] > votes["GRCh38"] else "GRCh38"
        loser = votes["GRCh38"] if winner == "GRCh37" else votes["GRCh37"]
        if loser == 0:
            evidence["method"] = "contig_lengths"
            return winner, evidence
        evidence["method"] = "contig_lengths_mixed"
        return "unknown", evidence
    if ref_hint:
        evidence["method"] = "reference_line"
        return ref_hint, evidence
    return "unknown", evidence


def read_vcf_header(vcf_path) -> str:
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    lines = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            lines.append(line.rstrip("\n"))
    return "\n".join(lines)


def detect_build(vcf_path) -> tuple:
    return detect_build_from_header(read_vcf_header(vcf_path))


def vcf_uses_chr_prefix(vcf_path) -> bool:
    header = read_vcf_header(vcf_path)
    m = _CONTIG.search(header)
    if m:
        fields = dict(kv.split("=", 1) for kv in m.group(1).split(",") if "=" in kv)
        return fields.get("ID", "").lower().startswith("chr")
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t", 1)[0].lower().startswith("chr")
    return False


def detect_build_from_markers(vcf_path, positions_37: list, chain_19to38) -> tuple:
    """Headerless fallback: count how many of the given GRCh37 (chrom, pos)
    marker sites carry a record in the VCF, versus their GRCh38 positions
    obtained through the hg19→hg38 chain. Needs a tabix index. Returns
    (build, evidence); 'unknown' when both counts are tiny or equal."""
    from cyvcf2 import VCF
    from pyliftover import LiftOver

    lo = LiftOver(str(chain_19to38))
    vcf = VCF(str(vcf_path))
    chr_prefix = vcf_uses_chr_prefix(vcf_path)
    hits = {"GRCh37": 0, "GRCh38": 0}
    for chrom, pos in positions_37:
        norm = normalize_chrom(str(chrom))
        for build, p in (("GRCh37", pos), ("GRCh38", None)):
            if build == "GRCh38":
                conv = lo.convert_coordinate(chain_chrom(norm), int(pos) - 1)
                if not conv or normalize_chrom(conv[0][0]) != norm:
                    continue
                p = conv[0][1] + 1
            region = f"{styled_chrom(norm, chr_prefix)}:{p}-{p}"
            try:
                if any(v.POS == p for v in vcf(region)):
                    hits[build] += 1
            except Exception:
                continue
    evidence = {"method": "marker_probe", "hits": hits, "n_markers": len(positions_37)}
    if max(hits.values()) >= 5 and hits["GRCh37"] != hits["GRCh38"]:
        return ("GRCh37" if hits["GRCh37"] > hits["GRCh38"] else "GRCh38"), evidence
    return "unknown", evidence


# ── chain files ──────────────────────────────────────────────────────────────

def chain_path(ref_dir, name: str) -> Path:
    return Path(ref_dir) / "liftover" / f"{name}.over.chain.gz"


def ensure_chain(ref_dir, name: str, download: bool = True) -> Path:
    path = chain_path(ref_dir, name)
    if path.exists() and path.stat().st_size > 0:
        return path
    if not download:
        raise FileNotFoundError(f"Chain file missing: {path}")
    url = CHAIN_URLS[name]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(path)
    return path


# ── FASTA access (hg19, UCSC 'chr' names, .fai index) ───────────────────────

class FastaIndex:
    """Minimal random access to a FASTA through its .fai (samtools-style)."""

    def __init__(self, fasta_path):
        self.path = Path(fasta_path)
        fai = Path(str(fasta_path) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"{fai} not found (run `samtools faidx {fasta_path}`)")
        self.index = {}
        with open(fai) as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fh = open(self.path, "rb")

    def has(self, name: str) -> bool:
        return name in self.index

    def fetch(self, name: str, start1: int, length: int) -> str:
        """Bases [start1, start1+length) in 1-based coordinates, upper-case."""
        total, offset, linebases, linewidth = self.index[name]
        start0 = start1 - 1
        if start0 < 0 or start0 + length > total:
            return ""
        out = []
        pos = start0
        remaining = length
        while remaining > 0:
            line_no, col = divmod(pos, linebases)
            self._fh.seek(offset + line_no * linewidth + col)
            take = min(remaining, linebases - col)
            out.append(self._fh.read(take).decode())
            pos += take
            remaining -= take
        return "".join(out).upper()

    def close(self):
        self._fh.close()


# ── liftover ─────────────────────────────────────────────────────────────────

def _flip_gt(gt: str) -> str:
    return "".join("1" if c == "0" else "0" if c == "1" else c for c in gt)


def _swap_sample(sample: str, fmt_keys: list) -> str:
    """REF/ALT swapped on a biallelic site: flip GT alleles, reverse AD and PL."""
    vals = sample.split(":")
    for i, key in enumerate(fmt_keys):
        if i >= len(vals):
            break
        if key == "GT":
            vals[i] = _flip_gt(vals[i])
        elif key == "AD":
            parts = vals[i].split(",")
            if len(parts) == 2:
                vals[i] = ",".join(reversed(parts))
        elif key in ("PL", "GL"):
            parts = vals[i].split(",")
            if len(parts) == 3:
                vals[i] = ",".join(reversed(parts))
    return ":".join(vals)


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def liftover_vcf(in_vcf, out_vcf, chain_file, fasta_path=None, bcftools: str = "bcftools",
                 progress=None, target_build: str = "GRCh37") -> dict:
    """VCF liftover through a UCSC chain (bgzipped, sorted, tabix-indexed output).
    Default direction is GRCh38 → GRCh37 (the pipeline's build); pass
    target_build="GRCh38" with the hg19ToHg38 chain for the reverse. Returns stats.

    Rules per record:
      * contig not in 1-22/X/Y/MT → dropped (contigs_skipped)
      * no chain mapping / more than one mapping → dropped (unmapped / ambiguous)
      * multi-base REF whose end does not map contiguously on the same
        chromosome and strand → dropped (unmapped)
      * minus-strand block: SNV alleles complemented; indels dropped
        (indel_minus_strand) — re-anchoring them safely needs more than a chain
      * FASTA given: hg19 REF must match; if it equals the single ALT instead,
        REF/ALT are swapped and GT/AD/PL adjusted (ref_alt_swapped); other
        mismatches dropped (ref_mismatch)
    Chromosome naming style ('chr' or not) of the input is preserved."""
    from pyliftover import LiftOver

    in_vcf, out_vcf = Path(in_vcf), Path(out_vcf)
    if target_build not in ("GRCh37", "GRCh38"):
        raise ValueError(f"target_build must be GRCh37 or GRCh38, got {target_build!r}")
    target_lengths = GRCH37_LENGTHS if target_build == "GRCh37" else GRCH38_LENGTHS
    source_build = "GRCh38" if target_build == "GRCh37" else "GRCh37"
    lo = LiftOver(str(chain_file))
    fasta = FastaIndex(fasta_path) if fasta_path and Path(fasta_path).exists() else None
    chr_prefix = vcf_uses_chr_prefix(in_vcf)

    stats = {k: 0 for k in ("n_in", "n_out", "contigs_skipped", "unmapped", "ambiguous",
                            "indel_minus_strand", "ref_mismatch", "ref_alt_swapped",
                            "strand_flipped")}
    stats["fasta_checked"] = fasta is not None

    tmp_dir = Path(tempfile.mkdtemp(prefix="bluegen_liftover_"))
    tmp_vcf = tmp_dir / "lifted_unsorted.vcf"
    opener = gzip.open if in_vcf.suffix == ".gz" else open
    with opener(in_vcf, "rt") as src, open(tmp_vcf, "w") as dst:
        for line in src:
            if line.startswith("##"):
                if line.startswith("##contig="):
                    m = _CONTIG.match(line.strip())
                    fields = dict(kv.split("=", 1) for kv in m.group(1).split(",") if "=" in kv)
                    norm = normalize_chrom(fields.get("ID", ""))
                    if norm in target_lengths:
                        fields["length"] = str(target_lengths[norm])
                        fields.pop("assembly", None)
                        line = "##contig=<" + ",".join(f"{k}={v}" for k, v in fields.items()) + ">\n"
                    elif norm == "MT":
                        pass
                    else:
                        continue  # non-standard contig: its records are dropped below
                elif line.startswith("##reference="):
                    line = (f"##reference={target_build} (lifted by BlueGen from {source_build}: "
                            + line[len("##reference="):])
                dst.write(line)
                continue
            if line.startswith("#CHROM"):
                dst.write(f"##bluegen_liftover=chain={Path(chain_file).name};"
                          f"ref_check={'hg19_fasta' if fasta else 'none'}\n")
                dst.write(line)
                continue

            stats["n_in"] += 1
            if progress and stats["n_in"] % 500000 == 0:
                progress(stats["n_in"])
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            norm = normalize_chrom(chrom)
            if norm not in STANDARD_CHROMS:
                stats["contigs_skipped"] += 1
                continue
            conv = lo.convert_coordinate(chain_chrom(norm), pos - 1)
            if not conv:
                stats["unmapped"] += 1
                continue
            if len(conv) > 1:
                stats["ambiguous"] += 1
                continue
            new_chrom_chain, new_pos0, strand, _score = conv[0]
            new_norm = normalize_chrom(new_chrom_chain)
            if new_norm not in STANDARD_CHROMS:
                stats["unmapped"] += 1
                continue
            alts = alt.split(",")
            is_indel = len(ref) != 1 or any(len(a) != 1 for a in alts if a not in (".", "*"))
            ref_len = len(ref)

            if strand == "-":
                if is_indel:
                    stats["indel_minus_strand"] += 1
                    continue
                stats["strand_flipped"] += 1
                new_pos = new_pos0 + 1
                ref = _revcomp(ref)
                alts = [_revcomp(a) if a not in (".", "*") else a for a in alts]
            else:
                new_pos = new_pos0 + 1
                if ref_len > 1:
                    conv_end = lo.convert_coordinate(chain_chrom(norm), pos - 1 + ref_len - 1)
                    if (len(conv_end) != 1 or normalize_chrom(conv_end[0][0]) != new_norm
                            or conv_end[0][2] != "+" or conv_end[0][1] != new_pos0 + ref_len - 1):
                        stats["unmapped"] += 1
                        continue

            if fasta is not None:
                fname = chain_chrom(new_norm)
                if fasta.has(fname):
                    genome_ref = fasta.fetch(fname, new_pos, len(ref))
                    if genome_ref != ref.upper():
                        if (not is_indel and len(alts) == 1 and alts[0] != "."
                                and genome_ref == alts[0].upper()):
                            ref, alts = alts[0], [ref]
                            stats["ref_alt_swapped"] += 1
                            if len(cols) > 9:
                                fmt_keys = cols[8].split(":")
                                cols[9:] = [_swap_sample(s, fmt_keys) for s in cols[9:]]
                        else:
                            stats["ref_mismatch"] += 1
                            continue

            cols[0] = styled_chrom(new_norm, chr_prefix)
            cols[1] = str(new_pos)
            cols[3] = ref
            cols[4] = ",".join(alts)
            dst.write("\t".join(cols) + "\n")
            stats["n_out"] += 1

    if fasta is not None:
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

    stats["dropped_total"] = stats["n_in"] - stats["n_out"]
    stats["output"] = str(out_vcf)
    stats["chain"] = str(chain_file)
    return stats


def write_build_record(path, record: dict) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(record, fh, indent=2)


__all__ = [
    "GRCH37_LENGTHS", "GRCH38_LENGTHS", "CHAIN_URLS", "normalize_chrom", "chain_chrom",
    "detect_build", "detect_build_from_header", "detect_build_from_markers", "read_vcf_header",
    "vcf_uses_chr_prefix", "ensure_chain", "chain_path", "FastaIndex", "liftover_vcf",
    "write_build_record",
]
