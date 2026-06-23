#!/usr/bin/env python3
"""
PRS Computation via PLINK --score — multi-trait, multi-sample.

Computes polygenic risk scores using the weighted sum method:
    PRS_i = Σ (β_j × G_ij)

Where:
    β_j  = GWAS effect size for SNP j (from curated SNP database)
    G_ij = genotype dosage for individual i at SNP j (0, 1, 2)

Input:
    --snp-db       CSV with columns: chrom, pos, effect_allele, weight, trait_category
    --bfile        PLINK binary dataset prefix (qc/qc_filtered)
    --output-dir   Output directory (default: prs/)
    --plink        Path to PLINK binary (default: plink)

Output:
    prs/prs_raw.csv   Columns: individual_id, trait, prs_raw, n_snps, n_snps_used

Usage:
    python prs_plink_score.py \\
        --snp-db data/snp_database_annotated.csv \\
        --bfile qc/qc_filtered \\
        --output-dir prs/ \\
        --plink tools/plink
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
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

        n_written = 0
        with open(score_file, "w") as fh:
            for _, row in trait_snps.iterrows():
                chrom = str(row.get("chrom", "")).replace("chr", "")
                pos = str(row.get("pos", "")).strip()
                allele = str(row.get("effect_allele", "")).strip()
                weight = row.get("weight", "1.0")
                vid = f"{chrom}:{pos}"
                if chrom and pos and allele and vid in bim_ids:
                    try:
                        w = float(weight)
                        fh.write(vid + chr(9) + allele + chr(9) + str(w) + chr(10))
                        n_written += 1
                    except ValueError:
                        continue

        if n_written == 0:
            score_file.unlink(missing_ok=True)
            continue

        out_prefix = out / f"tmp_{safe}"
        subprocess.run(
            [
                plink,
                "--bfile",
                bfile,
                "--score",
                str(score_file),
                "1",
                "2",
                "3",
                "--out",
                str(out_prefix),
                "--allow-extra-chr",
                "--threads",
                str(threads),
                "--memory",
                str(memory),
            ],
            capture_output=True,
            timeout=300,
        )

        profile_path = Path(str(out_prefix) + ".profile")
        if profile_path.exists():
            prof = pd.read_csv(profile_path, sep=r"[\t ]+", dtype={"IID": str})
            for _, prow in prof.iterrows():
                all_results.append(
                    {
                        "individual_id": str(prow["IID"]),
                        "trait": trait,
                        "prs_raw": float(prow.get("SCORE", prow.get("SCORESUM", 0))),
                        "n_snps": int(prow.get("CNT", 0)),
                        "n_snps_used": int(prow.get("CNT2", 0)),
                    }
                )
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


def main():
    parser = argparse.ArgumentParser(description="PRS Computation via PLINK --score")
    parser.add_argument("--snp-db", required=True, help="Path to SNP database CSV")
    parser.add_argument("--bfile", required=True, help="PLINK bfile prefix (qc/qc_filtered)")
    parser.add_argument("--output-dir", default="prs/", help="Output directory")
    parser.add_argument("--plink", default="plink", help="Path to PLINK binary")
    parser.add_argument("--threads", type=int, default=4, help="PLINK threads")
    parser.add_argument("--memory", type=int, default=8000, help="PLINK memory (MB)")
    args = parser.parse_args()

    compute_prs_plink_score(
        snp_db=args.snp_db,
        bfile=args.bfile,
        output_dir=args.output_dir,
        plink=args.plink,
        threads=args.threads,
        memory=args.memory,
    )


if __name__ == "__main__":
    main()
