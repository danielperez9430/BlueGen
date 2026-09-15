"""Stage F v3 — curated panel scored jointly with 1000G (RELEASE_PLAN 3.0.3).

Same synthetic-world technique as test_pgs_population_calibration.py: tiny
reference and user datasets built from VCFs with a real PLINK, a 6-SNP panel
across two traits, expectations computed by hand. Pins the bug this replaces
(user averaged over their ALT-carrying sites only → lactose intolerance
z = +5.4 on the real genome) and the array site policy.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PIPELINE = REPO_ROOT / "prs_research_pipeline"
sys.path.insert(0, str(PIPELINE))

from bluegen.scoring import compute_prs_joint  # noqa: E402
from bluegen.calibration import PopulationCalibrationV2  # noqa: E402


def find_plink():
    local = REPO_ROOT / "tools" / "plink"
    return str(local) if local.exists() else shutil.which("plink")


PLINK = find_plink()
pytestmark = pytest.mark.skipif(PLINK is None, reason="needs a PLINK 1.9 binary")

# (chrom, pos, ref, alt, effect_allele, weight, trait)
SNPS = [
    ("1", 100, "A", "G", "G", 0.50, "Trait A"),
    ("1", 200, "C", "T", "C", 0.25, "Trait A"),   # effect allele == REF
    ("1", 300, "G", "A", "A", 0.10, "Trait A"),
    ("2", 100, "T", "C", "C", 0.80, "Trait B"),
    ("2", 200, "A", "C", "A", 0.40, "Trait B"),   # effect allele == REF
    ("2", 300, "C", "G", "G", 0.30, "Trait B"),
]
POPS = {"EUR": 20, "EAS": 20}
FREQ = {"EUR": [0.10, 0.50, 0.30, 0.20, 0.60, 0.40], "EAS": [0.60, 0.20, 0.70, 0.50, 0.10, 0.80]}


def make_reference(rng):
    ids, pops, rows = [], [], []
    for pop, n in POPS.items():
        for k in range(n):
            ids.append(f"{pop}{k:02d}")
            pops.append(pop)
            rows.append([rng.binomial(2, FREQ[pop][j]) for j in range(len(SNPS))])
    return ids, pops, np.array(rows)


def write_vcf(path, sample_ids, alt_counts, keep_sites=None):
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids)]
    for j, (chrom, pos, ref, alt, *_rest) in enumerate(SNPS):
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
    _, _, ref, alt, eff, _, _ = SNPS[j]
    return alt_count if eff == alt else 2 - alt_count


def expected_average(alt_counts_row, sites):
    """PLINK default SCORE: sum(w × dosage) / (2 × loci)."""
    return sum(SNPS[j][5] * effect_dosage(alt_counts_row[j], j) for j in sites) / (2 * len(sites))


@pytest.fixture
def world(tmp_path):
    rng = np.random.RandomState(11)
    ids, pops, ref_counts = make_reference(rng)
    write_vcf(tmp_path / "ref.vcf", ids, ref_counts)
    vcf_to_bfile(tmp_path / "ref.vcf", tmp_path / "ref", keep_allele_order=True)
    (tmp_path / "panel.txt").write_text(
        "sample\tpop\tsuper_pop\n" + "".join(f"{i}\t{p}x\t{p}\n" for i, p in zip(ids, pops)))
    panel = pd.DataFrame([{
        "rsid": f"rs{j}", "gene": "G", "trait_category": t, "effect_allele": eff, "reference_allele": ref,
        "effect_direction": "+", "weight": w, "chrom": c, "pos": p} for j, (c, p, ref, alt, eff, w, t) in enumerate(SNPS)])
    panel.to_csv(tmp_path / "snp_db.csv", index=False)
    return {"dir": tmp_path, "ids": ids, "pops": pops, "ref_counts": ref_counts}


def run_joint(w, user_counts, keep_sites, site_policy):
    """keep_sites: which SNP indices appear in the user's VCF (WGS: only non-ref
    sites; array: every genotyped site, hom-ref included)."""
    write_vcf(w["dir"] / "user.vcf", ["S1"], user_counts, keep_sites=keep_sites)
    vcf_to_bfile(w["dir"] / "user.vcf", w["dir"] / "user", keep_allele_order=False)
    out = w["dir"] / "prs"
    df = compute_prs_joint(str(w["dir"] / "snp_db.csv"), str(w["dir"] / "user"), str(w["dir"] / "ref"),
                           str(w["dir"] / "panel.txt"), output_dir=str(out), plink=PLINK,
                           site_policy=site_policy, threads=1, memory=512, min_variants_per_chrom=1)
    ref = pd.read_csv(out / "prs_reference_raw.csv")
    return df, ref


def test_wgs_policy_fills_absent_panel_sites_as_hom_ref(world):
    w = world
    user = np.array([[1, 0, 0, 2, 0, 0]])           # ALT at sites 0 and 3 only
    df, ref = run_joint(w, user, keep_sites=[0, 3], site_policy="wgs")
    a = df[df.trait == "Trait A"].iloc[0]
    b = df[df.trait == "Trait B"].iloc[0]
    assert a["n_snps"] == 3 and a["n_snps_used"] == 3 and a["n_snps_panel"] == 3
    assert a["prs_raw"] == pytest.approx(expected_average(user[0], [0, 1, 2]), abs=1e-5)
    assert b["prs_raw"] == pytest.approx(expected_average(user[0], [3, 4, 5]), abs=1e-5)
    # The old user-only path would have averaged Trait A over site 0 alone:
    assert a["prs_raw"] != pytest.approx(expected_average(user[0], [0]), abs=1e-5)
    assert set(ref.individual_id) == set(w["ids"]) and ref.trait.nunique() == 2
    k = w["ids"].index("EUR03")
    ref_a = ref[(ref.individual_id == "EUR03") & (ref.trait == "Trait A")].iloc[0]
    assert ref_a["prs_raw"] == pytest.approx(expected_average(w["ref_counts"][k], [0, 1, 2]), abs=1e-5)


def test_array_policy_scores_only_genotyped_sites_for_everyone(world):
    """Array: sites 0,1,3,4 genotyped (site 1 hom-ref, present as 0/0); sites 2 and
    5 not on the array → excluded from the user AND the reference."""
    w = world
    user = np.array([[1, 0, 0, 2, 1, 0]])
    df, ref = run_joint(w, user, keep_sites=[0, 1, 3, 4], site_policy="array")
    a = df[df.trait == "Trait A"].iloc[0]
    b = df[df.trait == "Trait B"].iloc[0]
    assert (a["n_snps"], a["n_snps_used"], a["n_snps_panel"]) == (2, 2, 3)
    assert a["prs_raw"] == pytest.approx(expected_average(user[0], [0, 1]), abs=1e-5)
    assert b["prs_raw"] == pytest.approx(expected_average(user[0], [3, 4]), abs=1e-5)
    k = w["ids"].index("EAS05")
    ref_b = ref[(ref.individual_id == "EAS05") & (ref.trait == "Trait B")].iloc[0]
    assert ref_b["prs_raw"] == pytest.approx(expected_average(w["ref_counts"][k], [3, 4]), abs=1e-5)
    assert (df.site_policy == "array").all()


def test_reference_rows_feed_stage_h_and_a_reference_member_is_average(world):
    """A reference individual re-run as the user: joint scoring + Stage H's
    per-run distributions must give |z| small (it is a member of EUR)."""
    w = world
    k = w["ids"].index("EUR07")
    counts = w["ref_counts"][k:k + 1]
    nonref = [j for j in range(len(SNPS)) if counts[0, j] > 0]
    df, ref = run_joint(w, counts, keep_sites=nonref, site_policy="wgs")
    for trait in ("Trait A", "Trait B"):
        sites = [j for j, s in enumerate(SNPS) if s[6] == trait]
        row = df[df.trait == trait].iloc[0]
        assert row["prs_raw"] == pytest.approx(expected_average(counts[0], sites), abs=1e-5)
    cal = PopulationCalibrationV2()
    dists = cal.build_reference_distributions(str(w["dir"] / "prs" / "prs_reference_raw.csv"),
                                              str(w["dir"] / "panel.txt"), str(w["dir"] / "dists"))
    assert set(dists) == {"Trait A", "Trait B"} and set(dists["Trait A"]) == {"EUR", "EAS"}
    anc = w["dir"] / "anc.json"
    anc.write_text('{"assigned_population": "EUR"}')
    res = cal.calibrate_sample(str(w["dir"] / "prs" / "prs_raw.csv"), str(anc), str(w["dir"] / "cal"),
                               sample_id="EUR07", distributions=dists)
    assert len(res) == 2 and all(abs(r.z_score_population) < 2.5 for r in res)


def test_work_files_are_cleaned_up(world):
    w = world
    run_joint(w, np.array([[1, 0, 0, 1, 0, 0]]), keep_sites=[0, 3], site_policy="wgs")
    assert not (w["dir"] / "prs" / "_joint_work").exists()
