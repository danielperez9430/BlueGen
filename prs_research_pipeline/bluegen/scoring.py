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
