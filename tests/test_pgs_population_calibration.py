"""PGS Catalog population calibration (RELEASE_PLAN.md 3.0.1).

Drives scripts/utils/pgs_population_calibrate.py with a real PLINK binary
on tiny synthetic data built from VCFs, so the expected numbers can be
recomputed by hand in the test:

- the user and every reference sample must be summed over the SAME
  variants, with the user's absent sites treated as homozygous reference
  (the bug being pinned: user summed over their ALT-carrying sites only,
  reference over everything → z = -131);
- the z-score must be taken against the super-population from the ancestry
  model, not EUR by default;
- chromosomes absent from the user's data must be dropped from both sides
  and reported as reduced coverage.

Skipped when no PLINK 1.9 binary is available (CI's plain test job); runs
inside the Docker image, where plink is on PATH.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "prs_research_pipeline" / "scripts" / "utils" / "pgs_population_calibrate.py"


def find_plink():
    local = REPO_ROOT / "tools" / "plink"
    if local.exists():
        return str(local)
    return shutil.which("plink")


PLINK = find_plink()
pytestmark = pytest.mark.skipif(PLINK is None, reason="needs a PLINK 1.9 binary (tools/plink or PATH)")

# ── synthetic world ──────────────────────────────────────────────────────────
# 6 SNPs, REF/ALT fixed; effect allele is sometimes REF, sometimes ALT.
SNPS = [  # (chrom, pos, ref, alt, effect_allele, weight)
    ("1", 100, "A", "G", "G", 0.50),
    ("1", 200, "C", "T", "C", -0.25),
    ("1", 300, "G", "A", "A", 0.10),
    ("2", 100, "T", "C", "C", 0.80),
    ("2", 200, "A", "C", "A", -0.40),
    ("2", 300, "C", "G", "G", 0.30),
]
POPS = {"EUR": 20, "EAS": 20, "AFR": 12}
# ALT allele frequency per population per SNP (deliberately different)
FREQ = {
    "EUR": [0.10, 0.50, 0.30, 0.20, 0.60, 0.40],
    "EAS": [0.60, 0.20, 0.70, 0.50, 0.10, 0.80],
    "AFR": [0.30, 0.30, 0.50, 0.70, 0.40, 0.20],
}


def make_reference(rng):
    """Returns (sample_ids, super_pop per sample, alt-count matrix [n_samples, n_snps])."""
    ids, pops, rows = [], [], []
    for pop, n in POPS.items():
        for k in range(n):
            ids.append(f"{pop}{k:02d}")
            pops.append(pop)
            rows.append([rng.binomial(2, FREQ[pop][j]) for j in range(len(SNPS))])
    return ids, pops, np.array(rows)


def write_vcf(path, sample_ids, alt_counts, keep_sites=None):
    """alt_counts: [n_samples, n_snps]. keep_sites: indices to emit (None = all)."""
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids)]
    for j, (chrom, pos, ref, alt, _, _) in enumerate(SNPS):
        if keep_sites is not None and j not in keep_sites:
            continue
        calls = "\t".join(gt[int(alt_counts[s, j])] for s in range(len(sample_ids)))
        lines.append(f"{chrom}\t{pos}\t{chrom}:{pos}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{calls}")
    path.write_text("\n".join(lines) + "\n")


def vcf_to_bfile(vcf, prefix, keep_allele_order):
    args = [PLINK, "--vcf", str(vcf), "--make-bed", "--out", str(prefix), "--allow-extra-chr"]
    if keep_allele_order:
        args.append("--keep-allele-order")
    r = subprocess.run(args, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr


def effect_dosage(alt_count, j):
    _, _, ref, alt, eff, _ = SNPS[j]
    return alt_count if eff == alt else 2 - alt_count


def expected_score(alt_counts_row, sites):
    return sum(SNPS[j][5] * effect_dosage(alt_counts_row[j], j) for j in sites)


@pytest.fixture
def world(tmp_path):
    rng = np.random.RandomState(7)
    ids, pops, ref_counts = make_reference(rng)
    write_vcf(tmp_path / "ref.vcf", ids, ref_counts)
    vcf_to_bfile(tmp_path / "ref.vcf", tmp_path / "ref", keep_allele_order=True)
    (tmp_path / "panel.txt").write_text(
        "sample\tpop\tsuper_pop\n" + "".join(f"{i}\t{p}x\t{p}\n" for i, p in zip(ids, pops)))

    # The user: carries ALT at sites 0, 3 (het) and 5 (hom-alt); hom-ref elsewhere.
    # As in a real WGS→PLINK conversion, ONLY the non-ref sites are present.
    user_counts = np.array([[1, 0, 0, 1, 0, 2]])
    write_vcf(tmp_path / "user.vcf", ["S1"], user_counts, keep_sites=[0, 3, 5])
    vcf_to_bfile(tmp_path / "user.vcf", tmp_path / "user", keep_allele_order=False)

    pgs = tmp_path / "pgs" / "PGS000001"
    pgs.mkdir(parents=True)
    (pgs / "PGS000001_clean.score").write_text(
        "".join(f"{c}:{p}\t{eff}\t{w}\n" for c, p, _, _, eff, w in SNPS))
    (pgs / "PGS000001_hmPOS_GRCh37.txt").write_text(
        "#trait_reported=Synthetic trait\nrsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
        + "".join(f"rs{j}\t{c}\t{p}\t{eff}\t{ref}\t{w}\n" for j, (c, p, ref, alt, eff, w) in enumerate(SNPS)))
    return {"dir": tmp_path, "ids": ids, "pops": pops, "ref_counts": ref_counts, "user_counts": user_counts}


def run_calibrate(w, ancestry=None, extra=()):
    args = [sys.executable, str(SCRIPT),
            "--bfile", str(w["dir"] / "ref"), "--pop-panel", str(w["dir"] / "panel.txt"),
            "--pgs-dir", str(w["dir"] / "pgs"), "--user-bfile", str(w["dir"] / "user"),
            "--output-dir", str(w["dir"] / "out"), "--plink", PLINK,
            "--threads", "1", "--memory", "512", "--min-variants-per-chrom", "1", *extra]
    if ancestry is not None:
        aj = w["dir"] / "ANCESTRY_MODEL.json"
        aj.write_text(json.dumps({"assigned_population": ancestry, "method": "test"}))
        args += ["--ancestry-json", str(aj)]
    r = subprocess.run(args, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    csv = pd.read_csv(w["dir"] / "out" / "pgs_calibrated.csv")
    report = json.loads((w["dir"] / "out" / "pgs_calibration_report.json").read_text())
    return csv, report, r.stdout


def ref_pop_stats(w, pop, sites):
    scores = [expected_score(w["ref_counts"][s], sites)
              for s, p in enumerate(w["pops"]) if p == pop]
    return float(np.mean(scores)), float(np.std(scores, ddof=1))


def test_user_is_scored_on_the_same_variants_with_hom_ref_fill(world):
    csv, report, _ = run_calibrate(world, ancestry="EAS")
    assert len(csv) == 1
    row = csv.iloc[0]
    sites = list(range(len(SNPS)))
    expected_user = expected_score(world["user_counts"][0], sites)
    assert row["n_snps_matched"] == 6 and row["n_snps"] == 6
    assert row["coverage"] == pytest.approx(1.0)
    assert row["sample_score"] == pytest.approx(expected_user, abs=1e-4)
    assert row["individual_id"] == "S1"


def test_z_score_uses_the_inferred_population_not_eur(world):
    csv, report, _ = run_calibrate(world, ancestry="EAS")
    row = csv.iloc[0]
    sites = list(range(len(SNPS)))
    mean, std = ref_pop_stats(world, "EAS", sites)
    expected_user = expected_score(world["user_counts"][0], sites)
    assert row["reference_population"] == "EAS"
    assert row["ancestry_source"] == "inferred"
    assert row["ref_mean"] == pytest.approx(mean, abs=1e-4)
    assert row["ref_std"] == pytest.approx(std, abs=1e-4)
    assert row["z_score"] == pytest.approx((expected_user - mean) / std, abs=2e-3)
    assert report["methodology"]["reference_population"] == "EAS"
    # EUR numbers must differ (the populations were built with different frequencies)
    eur_mean, _ = ref_pop_stats(world, "EUR", sites)
    assert abs(eur_mean - mean) > 1e-3
    z_by_pop = json.loads(row["z_by_population"])
    assert set(z_by_pop) == {"EUR", "EAS", "AFR"}


def test_missing_ancestry_falls_back_to_eur_and_says_so(world):
    csv, report, out = run_calibrate(world, ancestry=None)
    row = csv.iloc[0]
    assert row["reference_population"] == "EUR"
    assert row["ancestry_source"] == "fallback"
    assert report["methodology"]["ancestry_source"] == "fallback"
    assert "falling back to EUR" in out


def test_unknown_population_label_falls_back(world):
    csv, _, out = run_calibrate(world, ancestry="MARTIAN")
    assert csv.iloc[0]["reference_population"] == "EUR"
    assert csv.iloc[0]["ancestry_source"] == "fallback"


def test_uncovered_chromosomes_are_dropped_from_both_sides(world):
    """User VCF with chr1 sites only: chr2 variants must not be hom-ref
    filled (they were never sequenced) — they leave the joint set for the
    reference too, and coverage drops to 3/6 → unreliable."""
    w = world
    user_counts = np.array([[1, 0, 2, 0, 0, 0]])
    write_vcf(w["dir"] / "user.vcf", ["S1"], user_counts, keep_sites=[0, 2])
    vcf_to_bfile(w["dir"] / "user.vcf", w["dir"] / "user", keep_allele_order=False)
    w["user_counts"] = user_counts
    csv, report, _ = run_calibrate(w, ancestry="EUR")
    row = csv.iloc[0]
    sites = [0, 1, 2]
    assert row["n_snps_matched"] == 3
    assert row["coverage"] == pytest.approx(0.5)
    assert not bool(row["reliable"])
    assert report["methodology"]["chromosomes_used"] == ["1"]
    mean, std = ref_pop_stats(w, "EUR", sites)
    assert row["ref_mean"] == pytest.approx(mean, abs=1e-4)
    assert row["sample_score"] == pytest.approx(expected_score(user_counts[0], sites), abs=1e-4)


def test_reference_members_score_near_their_own_population(world):
    """A reference individual re-scored as 'the user' must land inside its
    population's distribution (|z| small) — the property the real chr22
    3-sample 1000G test VCF is meant to satisfy end to end."""
    w = world
    k = w["ids"].index("EUR03")
    counts = w["ref_counts"][k:k + 1]
    nonref = [j for j in range(len(SNPS)) if counts[0, j] > 0]
    write_vcf(w["dir"] / "user.vcf", ["EUR03"], counts, keep_sites=nonref)  # same IID as a reference sample
    vcf_to_bfile(w["dir"] / "user.vcf", w["dir"] / "user", keep_allele_order=False)
    csv, _, _ = run_calibrate(w, ancestry="EUR")
    row = csv.iloc[0]
    assert row["individual_id"] == "EUR03"      # not merged into the reference sample of the same name
    assert abs(row["z_score"]) < 2.5
    assert row["sample_score"] == pytest.approx(expected_score(counts[0], range(len(SNPS))), abs=1e-4)


def test_dist_json_keeps_the_legacy_shape_for_portability_script(world):
    _, report, _ = run_calibrate(world, ancestry="EUR")
    dist = json.loads((world["dir"] / "out" / "ref_distributions" / "PGS000001_dist.json").read_text())
    assert {"EUR", "EAS", "AFR"} <= set(dist)
    assert {"n", "mean", "std", "median", "p5", "p95", "min", "max"} <= set(dist["EUR"])
    assert not (world["dir"] / "out" / "_joint.bed").exists(), "work files must be cleaned up"
