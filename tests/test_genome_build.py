"""Genome-build detection and GRCh38→GRCh37 liftover (RELEASE_PLAN 3.0.2).

Synthetic world: a tiny "hg19" FASTA (chr1, chr2; 1000 bp each), a UCSC-style
chain mapping an hg38 block on chr1 to hg19 on the + strand and an hg38
block on chr2 to hg19 on the − strand, and a hand-written hg38 VCF whose
records exercise every rule of liftover_vcf(). Expected outputs are computed
by hand from the chain geometry and the FASTA, not by running the code.

The round-trip test with the real UCSC chains and the chr22 3-sample VCF
runs only when both chain files are present under reference/liftover/
(never in CI).
"""

import gzip
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PIPELINE = REPO_ROOT / "prs_research_pipeline"
sys.path.insert(0, str(PIPELINE))

from bluegen import genome_build as gb  # noqa: E402

HAS_BCFTOOLS = shutil.which("bcftools") is not None
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


# ── detection ────────────────────────────────────────────────────────────────

def _header(contigs, ref_line=None):
    lines = ["##fileformat=VCFv4.2"]
    if ref_line:
        lines.append(ref_line)
    for cid, length in contigs:
        lines.append(f"##contig=<ID={cid},length={length}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(lines)


def test_normalize_chrom():
    assert gb.normalize_chrom("chr1") == "1"
    assert gb.normalize_chrom("1") == "1"
    assert gb.normalize_chrom("chrM") == "MT"
    assert gb.normalize_chrom("M") == "MT"
    assert gb.normalize_chrom("chrX") == "X"
    assert gb.normalize_chrom("chr6_ssto_hap7") == "6_ssto_hap7"


def test_detects_grch37_from_contig_lengths():
    build, ev = gb.detect_build_from_header(_header([(f"chr{c}", gb.GRCH37_LENGTHS[str(c)]) for c in (1, 2, 22)]))
    assert build == "GRCh37" and ev["method"] == "contig_lengths"


def test_detects_grch38_from_contig_lengths_without_chr_prefix():
    build, ev = gb.detect_build_from_header(_header([(str(c), gb.GRCH38_LENGTHS[str(c)]) for c in (1, 2, 3, 22)]))
    assert build == "GRCh38" and ev["votes"] == {"GRCh37": 0, "GRCh38": 4}


def test_mixed_contig_lengths_are_unknown():
    contigs = [("1", gb.GRCH37_LENGTHS["1"]), ("2", gb.GRCH38_LENGTHS["2"]), ("3", gb.GRCH37_LENGTHS["3"])]
    build, ev = gb.detect_build_from_header(_header(contigs))
    assert build == "unknown" and ev["method"] == "contig_lengths_mixed"


def test_reference_line_used_when_no_contig_lengths():
    build, ev = gb.detect_build_from_header(_header([], ref_line="##reference=file:///ref/hs37d5.fa"))
    assert build == "GRCh37" and ev["method"] == "reference_line"
    build, _ = gb.detect_build_from_header(_header([], ref_line="##reference=GRCh38_full_analysis_set.fa"))
    assert build == "GRCh38"


def test_no_signal_is_unknown():
    build, ev = gb.detect_build_from_header(_header([]))
    assert build == "unknown" and ev["method"] is None


def test_contig_lengths_beat_a_contradicting_reference_line():
    contigs = [(str(c), gb.GRCH38_LENGTHS[str(c)]) for c in (1, 2, 3)]
    build, _ = gb.detect_build_from_header(_header(contigs, ref_line="##reference=hg19"))
    assert build == "GRCh38"


@pytest.mark.skipif(not (PIPELINE / "test_samples_3eur.vcf.gz").exists(), reason="test VCF not present")
def test_real_test_vcf_is_grch37():
    build, ev = gb.detect_build(PIPELINE / "test_samples_3eur.vcf.gz")
    assert build == "GRCh37"


# ── synthetic liftover world ─────────────────────────────────────────────────

@pytest.fixture
def world(tmp_path):
    rng = np.random.RandomState(3)
    seqs = {"chr1": "".join(rng.choice(list("ACGT"), 1000)), "chr2": "".join(rng.choice(list("ACGT"), 1000))}
    fasta = tmp_path / "hg19.fa"
    fai_lines = []
    with open(fasta, "wb") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n".encode())
            offset = fh.tell()
            for i in range(0, len(seq), 50):
                fh.write((seq[i:i + 50] + "\n").encode())
            fai_lines.append(f"{name}\t{len(seq)}\t{offset}\t50\t51")
    (tmp_path / "hg19.fa.fai").write_text("\n".join(fai_lines) + "\n")

    # hg38 chr1 [100,200) → hg19 chr1 [300,400) on +   (offset +200)
    # hg38 chr2 [100,110) → hg19 chr2 reverse [500,510): forward pos0 = 999 - (500+off) = 499 - off
    chain = tmp_path / "hg38ToHg19.chain"
    chain.write_text(
        "chain 1000 chr1 1000 + 100 200 chr1 1000 + 300 400 1\n100\n\n"
        "chain 900 chr2 1000 + 100 110 chr2 1000 - 500 510 2\n10\n\n"
    )
    return {"dir": tmp_path, "seqs": seqs, "fasta": fasta, "chain": chain}


def hg19_base(world, chrom, pos1):
    return world["seqs"][chrom][pos1 - 1]


def other_base(b):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def test_fasta_index_fetch_crosses_line_boundaries(world):
    fa = gb.FastaIndex(world["fasta"])
    assert fa.fetch("chr1", 49, 4) == world["seqs"]["chr1"][48:52]
    assert fa.fetch("chr2", 1, 3) == world["seqs"]["chr2"][:3]
    assert fa.fetch("chr1", 999, 5) == ""  # past the end
    fa.close()


FORMAT_HEADERS = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">',
]


def _write_vcf(path, records, contig_lengths):
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh38_full_analysis_set_plus_decoy_hla.fa"] + FORMAT_HEADERS
    for cid, length in contig_lengths:
        lines.append(f"##contig=<ID={cid},length={length}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
    lines += ["\t".join(str(x) for x in r) for r in records]
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def _read_records(path):
    out = []
    with gzip.open(path, "rt") as fh:
        header = []
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            out.append(line.rstrip("\n").split("\t"))
    return header, out


@pytest.mark.skipif(not HAS_BCFTOOLS, reason="needs bcftools for sort/index")
def test_liftover_applies_every_rule(world):
    w = world
    b351 = hg19_base(w, "chr1", 351)          # hg38 chr1:151 → hg19 chr1:351 (+)
    b361 = hg19_base(w, "chr1", 361)          # swap case
    b371 = hg19_base(w, "chr1", 371)          # mismatch case
    del_ref = w["seqs"]["chr1"][380:382]      # hg19 chr1:381-382 (+), 2-base deletion REF
    b500 = hg19_base(w, "chr2", 500)          # hg38 chr2:101 → hg19 chr2:500 (−)
    b391 = hg19_base(w, "chr1", 391)          # multi-allelic case

    records = [
        ["chr1", 151, "a_ok", b351, other_base(b351), ".", "PASS", ".", "GT:AD:PL", "0/1:5,7:10,0,50"],
        ["chr1", 161, "b_swap", other_base(b361), b361, ".", "PASS", ".", "GT:AD:PL", "1/1:5,7:10,0,50"],
        ["chr1", 171, "c_mismatch", other_base(b371), other_base(other_base(b371)), ".", "PASS", ".", "GT", "0/1"],
        ["chr1", 181, "d_del", del_ref, del_ref[0], ".", "PASS", ".", "GT", "0/1"],
        ["chr2", 101, "e_minus", COMP[b500], COMP[other_base(b500)], ".", "PASS", ".", "GT", "0/1"],
        ["chr2", 105, "f_minus_indel", "AC", "A", ".", "PASS", ".", "GT", "0/1"],
        ["chr1", 50, "g_unmapped", "A", "G", ".", "PASS", ".", "GT", "0/1"],
        ["chr1_random", 10, "h_contig", "A", "G", ".", "PASS", ".", "GT", "0/1"],
        ["chr1", 191, "i_multi", other_base(b391), f"{b391},{other_base(other_base(b391))}", ".", "PASS", ".", "GT", "0/1"],
    ]
    src = w["dir"] / "in_hg38.vcf.gz"
    _write_vcf(src, records, [("chr1", gb.GRCH38_LENGTHS["1"]), ("chr2", gb.GRCH38_LENGTHS["2"]),
                              ("chr1_random", 12345)])
    out = w["dir"] / "out_hg19.vcf.gz"
    stats = gb.liftover_vcf(src, out, w["chain"], w["fasta"])

    assert stats["n_in"] == 9
    assert stats["n_out"] == 4
    assert stats["contigs_skipped"] == 1
    assert stats["unmapped"] == 1
    assert stats["indel_minus_strand"] == 1
    assert stats["ref_mismatch"] == 2          # c_mismatch + i_multi (no multi-allelic swap)
    assert stats["ref_alt_swapped"] == 1
    assert stats["strand_flipped"] == 1
    assert stats["fasta_checked"] is True
    assert (w["dir"] / "out_hg19.vcf.gz.tbi").exists()

    header, recs = _read_records(out)
    by_id = {r[2]: r for r in recs}
    assert [r[2] for r in recs] == ["a_ok", "b_swap", "d_del", "e_minus"]  # sorted by chrom/pos
    a = by_id["a_ok"]
    assert (a[0], a[1], a[3], a[4], a[9]) == ("chr1", "351", b351, other_base(b351), "0/1:5,7:10,0,50")
    b = by_id["b_swap"]
    assert (b[1], b[3], b[4]) == ("361", b361, other_base(b361))
    assert b[9] == "0/0:7,5:50,0,10"           # GT flipped, AD and PL reversed
    d = by_id["d_del"]
    assert (d[1], d[3], d[4]) == ("381", del_ref, del_ref[0])
    e = by_id["e_minus"]
    assert (e[0], e[1], e[3], e[4]) == ("chr2", "500", b500, other_base(b500))

    # Header rewritten to GRCh37: lengths, provenance, non-standard contig gone
    htext = "\n".join(header)
    assert f"##contig=<ID=chr1,length={gb.GRCH37_LENGTHS['1']}>" in htext
    assert "chr1_random" not in htext
    assert "##bluegen_liftover=chain=hg38ToHg19.chain;ref_check=hg19_fasta" in htext
    assert gb.detect_build(out)[0] == "GRCh37"


@pytest.mark.skipif(not HAS_BCFTOOLS, reason="needs bcftools for sort/index")
def test_liftover_without_fasta_keeps_alleles_and_reports_it(world):
    w = world
    b351 = hg19_base(w, "chr1", 351)
    records = [["chr1", 151, "x", other_base(b351), b351, ".", "PASS", ".", "GT", "1/1"]]
    src = w["dir"] / "in.vcf.gz"
    _write_vcf(src, records, [("chr1", gb.GRCH38_LENGTHS["1"])])
    out = w["dir"] / "out.vcf.gz"
    stats = gb.liftover_vcf(src, out, w["chain"], fasta_path=None)
    assert stats["fasta_checked"] is False and stats["n_out"] == 1 and stats["ref_alt_swapped"] == 0
    _, recs = _read_records(out)
    assert recs[0][3] == other_base(b351) and recs[0][9] == "1/1"   # untouched


@pytest.mark.skipif(not HAS_BCFTOOLS, reason="needs bcftools for index")
def test_marker_probe_detects_grch38_from_records(world):
    """Headerless VCF: six panel markers (GRCh37 positions) whose GRCh38
    positions carry records → GRCh38; the reverse chain is written on the fly."""
    import subprocess
    w = world
    chain_19to38 = w["dir"] / "hg19ToHg38.chain"
    chain_19to38.write_text("chain 1000 chr1 1000 + 300 400 chr1 1000 + 100 200 1\n100\n\n")
    markers = [("1", 351 + 5 * k) for k in range(6)]                  # hg19
    records = [["chr1", p - 200, f"m{k}", "A", "G", ".", "PASS", ".", "GT", "0/1"]  # hg38 = hg19 − 200
               for k, (_, p) in enumerate(markers)]
    src = w["dir"] / "headerless_hg38.vcf.gz"
    lines = ["##fileformat=VCFv4.2"] + FORMAT_HEADERS + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    lines += ["\t".join(str(x) for x in r) for r in records]
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(src)], input="\n".join(lines) + "\n",
                   text=True, check=True, capture_output=True)
    subprocess.run(["bcftools", "index", "-t", str(src)], check=True, capture_output=True)
    assert gb.detect_build(src)[0] == "unknown"
    build, ev = gb.detect_build_from_markers(src, markers, chain_19to38)
    assert build == "GRCh38" and ev["hits"] == {"GRCh37": 0, "GRCh38": 6}


# ── real chains: round trip on the chr22 test VCF (local only) ──────────────

REAL_CHAINS = [gb.chain_path(PIPELINE / "reference", n) for n in ("hg19ToHg38", "hg38ToHg19")]
TEST_VCF = PIPELINE / "test_samples_3eur.vcf.gz"
HG19_FA = PIPELINE / "reference" / "hg19" / "hg19.fa"


@pytest.mark.skipif(not (HAS_BCFTOOLS and all(p.exists() for p in REAL_CHAINS) and TEST_VCF.exists()),
                    reason="needs the UCSC chain files under reference/liftover/ and the test VCF")
def test_round_trip_chr22_with_real_chains(tmp_path):
    """hg19 → hg38 (no FASTA) → hg19 (hg19 FASTA when available). Nearly every
    record must come back at its original position with original alleles."""
    to38 = tmp_path / "hg38.vcf.gz"
    back = tmp_path / "back_hg19.vcf.gz"
    s1 = gb.liftover_vcf(TEST_VCF, to38, REAL_CHAINS[0], fasta_path=None, target_build="GRCh38")
    assert s1["n_out"] > 0.98 * s1["n_in"]
    s2 = gb.liftover_vcf(to38, back, REAL_CHAINS[1], HG19_FA if HG19_FA.exists() else None)
    assert s2["n_out"] > 0.98 * s2["n_in"]
    _, orig = _read_records(TEST_VCF)
    _, rt = _read_records(back)
    orig_keys = {(r[0], r[1], r[3], r[4]) for r in orig}
    rt_keys = {(r[0], r[1], r[3], r[4]) for r in rt}
    assert len(rt_keys & orig_keys) / len(orig_keys) > 0.97
    assert gb.detect_build(to38)[0] == "GRCh38"
    assert gb.detect_build(back)[0] == "GRCh37"
    assert json.dumps(s2)  # stats are JSON-serialisable for input_build.json
