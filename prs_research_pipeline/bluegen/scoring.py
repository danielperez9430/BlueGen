"""
Stage F — PRS Computation via PLINK --score (IMPROVEMENT_PLAN.md 2.1).

Extracted from scripts/prs/prs_plink_score.py, which is now a thin CLI
wrapper around this module. Logic unchanged.

Computes polygenic risk scores using the weighted sum method:
    PRS_i = Sum(beta_j * G_ij)

Where:
    beta_j = GWAS effect size for SNP j (from curated SNP database)
    G_ij   = genotype dosage for individual i at SNP j (0, 1, 2)

This is the only one of the three critical-path stages that shells out to a
real external binary (PLINK). The SNP-matching logic (which panel rows match
a given set of .bim variant ids) is split out as build_score_rows() - a pure
function with no I/O and no subprocess, fully unit-testable with a fake
bim_ids set. The PLINK subprocess call and .profile parsing stay as
separate, narrowly-scoped functions so they can be tested independently
(e.g. with a mocked subprocess.run) without needing a real PLINK binary.
"""

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def fix_duplicate_bim_ids(bim_path: Path, output_dir: Path, bfile: str) -> str:
    """Fix duplicate variant IDs in BIM file. Returns (possibly new) bfile prefix."""
    bim = pd.read_csv(bim_path, sep=r"\s+", header=None, dtype=str)
    bim.columns = ["chr", "vid", "cm", "pos", "a1", "a2"]
    dups = bim[bim["vid"].duplicated()]["vid"].unique()

    if len(dups) == 0:
        return bfile

    dup_counts = {}
    new_vids = []
    for _, row in bim.iterrows():
        v = row["vid"]
        if v in dups:
            if v not in dup_counts:
                dup_counts[v] = 1
                new_vids.append(v)
            else:
                dup_counts[v] += 1
                new_vids.append(f"{v}_{dup_counts[v]}")
        else:
            new_vids.append(v)

    bim["vid"] = new_vids
    tmp_bim = output_dir / "qc_dedup.bim"
    bim.to_csv(tmp_bim, sep="\t", header=False, index=False)
    for ext in [".bed", ".fam"]:
        shutil.copy2(Path(bfile + ext), output_dir / ("qc_dedup" + ext))

    print(f"  PRS: fixed {len(dups)} duplicate bim IDs")
    return str(output_dir / "qc_dedup")


def build_score_rows(trait_snps: pd.DataFrame, bim_ids: set) -> List[Tuple[str, str, float]]:
    """
    Pure SNP-matching logic for one trait, no I/O, no PLINK.

    Matches panel rows to genotyped variants by chrom:pos (same key PLINK's
    .bim uses), same as the rest of this pipeline - a SNP whose position
    isn't in bim_ids is silently dropped (documented existing behavior,
    see tests/test_snp_positions.py / test_allele_strand_consistency.py,
    not changed here).

    Args:
        trait_snps: rows of the curated SNP panel for one trait_category,
            with at least chrom/pos/effect_allele/weight columns (str dtype,
            matching pd.read_csv(snp_db, dtype=str) upstream), plus an
            optional effect_direction column ('+'/'-', defaults to '+').
        bim_ids: set of "chrom:pos" variant ids present in the genotype data.

    Returns:
        List of (vid, effect_allele, signed_weight) tuples ready to write to
        a PLINK --score file, for SNPs that matched and had a parseable
        weight. PLINK's --score format has no separate direction column, so
        effect_direction='-' (a protective/risk-lowering effect_allele) is
        applied here by negating the weight - this was previously silently
        ignored (IMPROVEMENT_PLAN.md follow-up, found via a full-panel
        polarity audit): every row used effect_direction='+' by construction
        except 3, which were being added to the score instead of subtracted.
    """
    rows = []
    for _, row in trait_snps.iterrows():
        chrom = str(row.get("chrom", "")).replace("chr", "")
        pos = str(row.get("pos", "")).strip()
        allele = str(row.get("effect_allele", "")).strip()
        weight = row.get("weight", "1.0")
        direction = str(row.get("effect_direction", "+")).strip()
        vid = f"{chrom}:{pos}"
        if chrom and pos and allele and vid in bim_ids:
            try:
                w = float(weight)
                if direction == "-":
                    w = -w
                rows.append((vid, allele, w))
            except ValueError:
                continue
    return rows


def write_score_file(score_path: Path, rows: List[Tuple[str, str, float]]) -> None:
    """Write PLINK --score input: tab-separated vid, effect_allele, weight."""
    with open(score_path, "w") as fh:
        for vid, allele, weight in rows:
            fh.write(vid + chr(9) + allele + chr(9) + str(weight) + chr(10))


def run_plink_score(
    plink: str, bfile: str, score_file: Path, out_prefix: Path,
    threads: int = 4, memory: int = 8000,
) -> subprocess.CompletedProcess:
    """Invoke PLINK --score. Thin subprocess wrapper, mockable in tests."""
    return subprocess.run(
        [
            plink, "--bfile", bfile, "--score", str(score_file), "1", "2", "3",
            "--out", str(out_prefix), "--allow-extra-chr",
            "--threads", str(threads), "--memory", str(memory),
        ],
        capture_output=True, timeout=300,
    )


def parse_plink_profile(profile_path: Path, trait: str) -> List[Dict]:
    """Parse a PLINK .profile file into per-individual result dicts."""
    results = []
    prof = pd.read_csv(profile_path, sep=r"[\t ]+", dtype={"IID": str})
    for _, prow in prof.iterrows():
        results.append({
            "individual_id": str(prow["IID"]),
            "trait": trait,
            "prs_raw": float(prow.get("SCORE", prow.get("SCORESUM", 0))),
            "n_snps": int(prow.get("CNT", 0)),
            "n_snps_used": int(prow.get("CNT2", 0)),
        })
    return results


def compute_prs_plink_score(
    snp_db: str,
    bfile: str,
    output_dir: str = "prs/",
    plink: str = "plink",
    threads: int = 4,
    memory: int = 8000,
) -> pd.DataFrame:
    """
    Compute PRS per trait using PLINK --score.

    Returns DataFrame with columns: individual_id, trait, prs_raw, n_snps, n_snps_used.
    """
    db = pd.read_csv(snp_db, dtype=str)
    trait_col = "trait_category"
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    bim_path = Path(bfile + ".bim")
    bfile = fix_duplicate_bim_ids(bim_path, out, bfile)

    bim = pd.read_csv(Path(bfile + ".bim"), sep=r"\s+", header=None, dtype=str)
    bim.columns = ["chr", "vid", "cm", "pos", "a1", "a2"]
    bim_ids = set(bim["vid"].values)

    traits = db[trait_col].dropna().unique()
    all_results = []

    for trait in traits:
        trait_snps = db[db[trait_col] == trait]
        safe = (
            trait.lower()
            .replace(" ", "_")
            .replace("&", "and")
            .replace("/", "_")
            .replace("(", "")
            .replace(")", "")
        )
        score_file = out / f"tmp_{safe}.score"

        rows = build_score_rows(trait_snps, bim_ids)
        if not rows:
            continue
        write_score_file(score_file, rows)

        out_prefix = out / f"tmp_{safe}"
        run_plink_score(plink, bfile, score_file, out_prefix, threads=threads, memory=memory)

        profile_path = Path(str(out_prefix) + ".profile")
        if profile_path.exists():
            all_results.extend(parse_plink_profile(profile_path, trait))
            profile_path.unlink()

        score_file.unlink(missing_ok=True)
        for ext in [".log", ".nosex", ".nopred"]:
            Path(str(out_prefix) + ext).unlink(missing_ok=True)

    if all_results:
        df = pd.DataFrame(all_results)
        n_samples = df["individual_id"].nunique()
        n_traits = df["trait"].nunique()
        df.to_csv(out / "prs_raw.csv", index=False)
        print(
            f"  PRS: {n_traits} traits, {n_samples} samples, "
            f'{df["n_snps_used"].sum()}/{df["n_snps"].sum()} SNPs'
        )
    else:
        print("  PRS: 0 scores — no variants matched")
        df = pd.DataFrame()
        df.to_csv(out / "prs_raw.csv", index=False)

    for ext in [".bed", ".bim", ".fam"]:
        Path(str(out / "qc_dedup") + ext).unlink(missing_ok=True)

    return df


def compute_prs_joint(
    snp_db: str,
    user_bfile: str,
    ref_bfile: str,
    pop_panel: str,
    output_dir: str = "prs/",
    plink: str = "plink",
    site_policy: str = "wgs",
    threads: int = 4,
    memory: int = 8000,
    min_variants_per_chrom: int = 100,
) -> pd.DataFrame:
    """
    Stage F v3 (RELEASE_PLAN 3.0.3): score the user JOINTLY with the 1000 Genomes
    reference on one identical panel-variant set.

    compute_prs_plink_score() scored the user alone on the variants present in
    their own PLINK file. A WGS-derived file only holds sites where the sample
    carries an ALT allele, so a trait with 5 panel SNPs where the user is
    hom-ref at 3 was averaged over the 2 ALT-carrying sites while the
    reference distribution (Stage H) averaged all 5 → e.g. lactose intolerance
    z = +5.4, percentile 100. Here the user is merged into the reference
    (bluegen.joint) with hom-ref filling (site_policy="wgs") or restricted to
    genotyped sites (site_policy="array"), and every trait is scored with one
    PLINK --score for user and reference together (default per-allele average,
    the same units the reference distributions and PC betas always used).

    Writes:
      prs/prs_raw.csv            user rows: individual_id, trait, prs_raw, n_snps,
                                 n_snps_used, n_snps_panel, site_policy
      prs/prs_reference_raw.csv  reference rows: individual_id, trait, prs_raw —
                                 input for Stage H's per-run distributions and
                                 Stage G's per-run PC betas
    Returns the user DataFrame.
    """
    from .joint import (JointDatasetError, build_joint_dataset, remove_joint_dataset, read_bim_ids,
                        user_covered_chromosomes)

    db = pd.read_csv(snp_db, dtype=str)
    trait_col = "trait_category"
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_joint_work"

    panel_ids = set()
    for _, row in db.iterrows():
        chrom = str(row.get("chrom", "")).replace("chr", "").strip()
        pos = str(row.get("pos", "")).strip()
        if chrom and pos and chrom != "nan" and pos != "nan":
            panel_ids.add(f"{chrom}:{pos}")

    chroms = None
    if site_policy == "wgs":
        chroms = user_covered_chromosomes(str(user_bfile) + ".bim", min_variants_per_chrom)
        if not chroms:
            print("  PRS: user data covers no autosome with enough variants — scoring nothing")
            pd.DataFrame().to_csv(out / "prs_raw.csv", index=False)
            return pd.DataFrame()

    try:
        joint, user_map, info = build_joint_dataset(
            plink, ref_bfile, user_bfile, panel_ids, work, threads=str(threads), memory=str(memory),
            site_policy=site_policy, chroms=chroms, prefix="_panel_joint")
    except JointDatasetError as e:
        # e.g. an array with no probe on any panel SNP: no scores, but the
        # rest of the pipeline (ClinVar, PharmGKB, ancestry) still runs.
        print(f"  PRS: 0 scores — {e}")
        pd.DataFrame().to_csv(out / "prs_raw.csv", index=False)
        pd.DataFrame().to_csv(out / "prs_reference_raw.csv", index=False)
        return pd.DataFrame()
    joint_ids = read_bim_ids(str(joint) + ".bim")
    user_genotyped = info["user_bim_ids"]
    print(f"  PRS (joint, {site_policy}): {info['n_joint_variants']}/{len(panel_ids)} panel variants in the "
          f"joint set, {info['n_samples']} samples ({len(user_map)} user)"
          + (f", {info['n_merge_conflicts_excluded']} merge conflicts excluded"
             if info["n_merge_conflicts_excluded"] else ""))

    user_rows, ref_rows = [], []
    for trait in db[trait_col].dropna().unique():
        trait_snps = db[db[trait_col] == trait]
        rows = build_score_rows(trait_snps, joint_ids)
        if not rows:
            continue
        safe = "".join(c if c.isalnum() else "_" for c in str(trait).lower())
        score_file = work / f"tmp_{safe}.score"
        write_score_file(score_file, rows)
        out_prefix = work / f"tmp_{safe}"
        run_plink_score(plink, str(joint), score_file, out_prefix, threads=threads, memory=memory)
        profile_path = Path(str(out_prefix) + ".profile")
        if not profile_path.exists():
            continue
        prof = pd.read_csv(profile_path, sep=r"[\t ]+", dtype={"FID": str, "IID": str}, engine="python")
        n_genotyped = sum(1 for vid, _, _ in rows if vid in user_genotyped)
        for _, prow in prof.iterrows():
            iid = str(prow["IID"])
            score = float(prow.get("SCORE", prow.get("SCORESUM", 0)))
            n_loci = int(prow.get("CNT", 0)) // 2
            if iid in user_map:
                user_rows.append({
                    "individual_id": user_map[iid], "trait": trait, "prs_raw": score,
                    "n_snps": n_loci,
                    "n_snps_used": n_loci if site_policy == "wgs" else min(n_genotyped, n_loci),
                    "n_snps_panel": int(len(trait_snps)),
                    "site_policy": site_policy,
                })
            else:
                ref_rows.append({"individual_id": iid, "trait": trait, "prs_raw": score})
        profile_path.unlink(missing_ok=True)
        score_file.unlink(missing_ok=True)
        for ext in (".log", ".nosex", ".nopred"):
            Path(str(out_prefix) + ext).unlink(missing_ok=True)

    remove_joint_dataset(joint, extra_prefixes=(work / "_panel_joint_ref",))
    try:
        work.rmdir()
    except OSError:
        pass

    df = pd.DataFrame(user_rows)
    df.to_csv(out / "prs_raw.csv", index=False)
    pd.DataFrame(ref_rows).to_csv(out / "prs_reference_raw.csv", index=False)
    if user_rows:
        print(f"  PRS: {df['trait'].nunique()} traits, {df['individual_id'].nunique()} sample(s), "
              f"{int(df['n_snps_used'].sum())}/{int(df['n_snps'].sum())} SNPs genotyped "
              f"({int(df['n_snps_panel'].sum())} in panel); reference rows: {len(ref_rows)}")
    else:
        print("  PRS: 0 scores — no panel variant matched the joint set")
    return df
