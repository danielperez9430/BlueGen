"""Stage B QC on single-sample datasets (RELEASE_PLAN 3.0.3).

MAF and HWE are cohort statistics. On one sample every homozygous site has
MAF = 0, so --maf 0.01 deleted every hom-ALT variant of a WGS run (1.66 M of
4.12 M on the real genome) and every hom-REF call of an array file. Those
sites were then treated as absent (== hom-REF) by the joint scoring. The
stage now skips both filters below 50 samples unless --force-freq-filters.
Drives the real bash stage with a real PLINK on a 3-site single-sample VCF.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
STAGE_B = REPO_ROOT / "prs_research_pipeline" / "scripts" / "stages" / "02_quality_control.sh"


def find_plink():
    local = REPO_ROOT / "tools" / "plink"
    return str(local) if local.exists() else shutil.which("plink")


PLINK = find_plink()
pytestmark = pytest.mark.skipif(
    PLINK is None or sys.platform.startswith("win") or shutil.which("bash") is None,
    reason="needs PLINK 1.9 and bash",
)


def _single_sample_bfile(tmp_path):
    vcf = tmp_path / "one.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t100\thomref\tA\tG\t.\tPASS\t.\tGT\t0/0\n"
        "1\t200\thet\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
        "1\t300\thomalt\tG\tA\t.\tPASS\t.\tGT\t1/1\n"
    )
    r = subprocess.run([PLINK, "--vcf", str(vcf), "--make-bed", "--out", str(tmp_path / "cohort"),
                        "--allow-extra-chr", "--keep-allele-order"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return tmp_path / "cohort"


def _run_stage_b(bfile, out_dir, *extra):
    r = subprocess.run(["bash", str(STAGE_B), "--bfile", str(bfile), "--out-dir", str(out_dir),
                        "--plink", PLINK, "--threads", "1", "--memory", "256", *extra],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-500:]
    bim = out_dir / "qc_filtered.bim"
    return r.stdout, [line.split()[1] for line in bim.read_text().splitlines()]


def test_single_sample_keeps_homozygous_sites(tmp_path):
    bfile = _single_sample_bfile(tmp_path)
    out, ids = _run_stage_b(bfile, tmp_path / "qc")
    assert "MAF/HWE filters SKIPPED" in out
    assert set(ids) == {"homref", "het", "homalt"}


def test_force_freq_filters_reproduces_the_old_behaviour(tmp_path):
    bfile = _single_sample_bfile(tmp_path)
    out, ids = _run_stage_b(bfile, tmp_path / "qc", "--force-freq-filters")
    assert "MAF filter" in out
    assert ids == ["het"]           # exactly the pre-3.0.3 loss: only heterozygous sites survive


def test_threshold_is_configurable(tmp_path):
    bfile = _single_sample_bfile(tmp_path)
    out, ids = _run_stage_b(bfile, tmp_path / "qc", "--min-samples-for-freq-filters", "1")
    assert "MAF filter" in out and ids == ["het"]
