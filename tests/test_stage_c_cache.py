"""Stage C (03_ld_ancestry_prune.sh) prune-set cache (RELEASE_PLAN.md 2.1.6).

The ancestry-matched prune set is derived only from the 1000G reference and
the --indep-pairwise / --union-mode parameters, so it is cached next to the
reference and reused across samples (45 min -> 30 s). Before 2.1.0 the cache
file name carried no parameters, so a run with e.g. --r2 0.5 silently reused
a set pruned at r2=0.2. These tests drive the real bash script with a stub
`plink` on PATH: the cache-hit path only needs `plink --extract`, and the
cache-miss path is aborted at its first PLINK call (the 1000G QC step) so we
can assert the cache was *not* used without needing real reference data.
"""

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "prs_research_pipeline" / "scripts" / "stages" / "03_ld_ancestry_prune.sh"

pytestmark = pytest.mark.skipif(
    sys.platform.startswith("win") or shutil.which("bash") is None,
    reason="needs bash",
)

STUB_PLINK = """#!/usr/bin/env bash
# Minimal PLINK stand-in for the cache tests.
out=""; extract=""; qc=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) out="$2"; shift 2 ;;
        --extract) extract="$2"; shift 2 ;;
        --geno) qc=1; shift 2 ;;
        *) shift ;;
    esac
done
if [[ $qc == 1 ]]; then
    # First call of the cache-MISS path (1000G QC). Leave a marker and fail
    # so the script stops here instead of needing pandas + real 1000G data.
    touch "$(dirname "$out")/QC_CALLED"
    exit 1
fi
if [[ -n "$extract" ]]; then
    cp "$extract" "${out}.bim"
    : > "${out}.bed"; : > "${out}.fam"
fi
exit 0
"""


@pytest.fixture
def sandbox(tmp_path):
    ref = tmp_path / "ref" / "1000G_full"
    ref.parent.mkdir()
    for ext in ("bed", "bim", "fam"):
        (ref.parent / f"1000G_full.{ext}").write_text("")
    target = tmp_path / "qc_filtered"
    (tmp_path / "qc_filtered.bed").write_text("")
    (tmp_path / "qc_filtered.bim").write_text(
        "\n".join(f"1\trs{i}\t0\t{i * 1000}\tA\tG" for i in range(1, 11)) + "\n"
    )
    (tmp_path / "panel.txt").write_text("sample\tpop\tsuper_pop\n")
    stub = tmp_path / "bin" / "plink"
    stub.parent.mkdir()
    stub.write_text(STUB_PLINK)
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
    return {"root": tmp_path, "ref": ref, "target": target, "plink": stub}


def run_stage_c(sb, *extra):
    cmd = [
        "bash", str(SCRIPT),
        "--bfile", str(sb["target"]),
        "--1000g-bfile", str(sb["ref"]),
        "--population-panel", str(sb["root"] / "panel.txt"),
        "--out-dir", str(sb["root"] / "out"),
        "--plink", str(sb["plink"]),
        *extra,
    ]
    return subprocess.run(cmd, cwd=sb["root"], capture_output=True, text=True,
                          env={**os.environ, "LC_ALL": "C"}, timeout=60)


def _write_prune_set(path, snps):
    path.write_text("\n".join(snps) + "\n")


def test_cache_hit_with_matching_parameters_is_used(sandbox):
    cache = Path(f"{sandbox['ref']}_ancestry_pruned_w50_s5_r0.2_conservative.txt")
    _write_prune_set(cache, ["rs1", "rs3", "rs5"])
    r = run_stage_c(sandbox)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Using cached ancestry-pruned SNP list" in r.stdout
    pruned = (sandbox["root"] / "out" / "ld_pruned_dataset.bim").read_text().split()
    assert pruned == ["rs1", "rs3", "rs5"]
    assert not (sandbox["root"] / "out" / "1000G_ld_reference" / "QC_CALLED").exists()


def test_legacy_cache_is_migrated_and_used_for_default_parameters(sandbox):
    legacy = Path(f"{sandbox['ref']}_ancestry_pruned_snps.txt")
    _write_prune_set(legacy, ["rs2", "rs4"])
    r = run_stage_c(sandbox)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Migrating legacy prune cache" in r.stdout
    assert Path(f"{sandbox['ref']}_ancestry_pruned_w50_s5_r0.2_conservative.txt").exists()
    assert "Using cached ancestry-pruned SNP list" in r.stdout
    assert legacy.exists(), "migration must copy, not move, the legacy file"


def test_cache_with_different_parameters_is_not_reused(sandbox):
    """A default-params cache (legacy or new) must NOT satisfy a --r2 0.5 run."""
    _write_prune_set(Path(f"{sandbox['ref']}_ancestry_pruned_snps.txt"), ["rs2"])
    _write_prune_set(Path(f"{sandbox['ref']}_ancestry_pruned_w50_s5_r0.2_conservative.txt"), ["rs2"])
    r = run_stage_c(sandbox, "--r2", "0.5")
    assert r.returncode != 0  # stub aborts at the QC step of the miss path
    assert "Using cached ancestry-pruned SNP list" not in r.stdout
    assert "Migrating legacy prune cache" not in r.stdout
    assert (sandbox["root"] / "out" / "1000G_ld_reference" / "QC_CALLED").exists()
    assert not Path(f"{sandbox['ref']}_ancestry_pruned_w50_s5_r0.5_conservative.txt").exists()


def test_union_mode_is_part_of_the_cache_key(sandbox):
    _write_prune_set(Path(f"{sandbox['ref']}_ancestry_pruned_w50_s5_r0.2_inclusive.txt"), ["rs9"])
    r = run_stage_c(sandbox, "--union-mode", "inclusive")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Using cached ancestry-pruned SNP list" in r.stdout
    pruned = (sandbox["root"] / "out" / "ld_pruned_dataset.bim").read_text().split()
    assert pruned == ["rs9"]
