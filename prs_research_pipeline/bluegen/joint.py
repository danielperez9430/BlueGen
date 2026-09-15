"""
Joint user + 1000 Genomes PLINK dataset builder (RELEASE_PLAN 3.0.1 / 3.0.3).

Both polygenic-score paths (the curated panel in Stage F and the PGS Catalog
calibration) must score the user and every reference individual over the
SAME variants with the SAME allele orientation, otherwise their sums are not
comparable. This module builds that dataset once per set of needed variants:

  1. extract the needed variant IDs (chrom:pos) from the reference and from
     the user, optionally restricted to autosomes the user's data covers;
  2. rename the user samples (USER_<iid>) so a 1000G member re-run as the
     user can never be merged into its own reference record;
  3. merge, reference first, with --keep-allele-order so A2 stays the hg19
     REF allele of the 1000G conversion (.missnp conflicts are excluded from
     both sides and the merge retried);
  4. site policy:
       wgs   — --fill-missing-a2: a site absent from a WGS-derived VCF is a
               homozygous-reference call, so it is filled as such (the same
               assumption 06_prs_compute.py makes);
       array — the needed set is first intersected with the sites the user
               was actually genotyped on; nothing is filled, the few
               genuinely missing calls stay missing (PLINK mean-imputes them
               inside --score).

Returns the joint prefix, the {USER_iid: iid} map and a small info dict.
"""

from __future__ import annotations

import sys
from pathlib import Path
import subprocess

AUTOSOMES = [str(c) for c in range(1, 23)]
SITE_POLICIES = ("wgs", "array")


def run_plink(plink, args, timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run([str(plink)] + [str(a) for a in args], capture_output=True, text=True, timeout=timeout)


def read_bim_ids(bim_path) -> set:
    ids = set()
    with open(bim_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                ids.add(parts[1])
    return ids


def read_fam_ids(fam_path) -> list:
    ids = []
    with open(fam_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                ids.append((parts[0], parts[1]))
    return ids


def user_covered_chromosomes(user_bim, min_variants: int) -> list:
    """Autosomes on which the user's dataset has at least `min_variants`
    variants — hom-ref filling is only valid where the sample was called."""
    counts = {}
    with open(user_bim) as fh:
        for line in fh:
            chrom = line.split("\t", 1)[0].split(" ", 1)[0]
            counts[chrom] = counts.get(chrom, 0) + 1
    return [c for c in AUTOSOMES if counts.get(c, 0) >= min_variants]


def rename_user_samples(fam: Path) -> dict:
    rows, mapping = [], {}
    with open(fam) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6:
                continue
            new_iid = f"USER_{parts[1]}"
            mapping[new_iid] = parts[1]
            rows.append("\t".join(["USER", new_iid] + parts[2:]))
    fam.write_text("\n".join(rows) + "\n")
    return mapping


class JointDatasetError(RuntimeError):
    """The joint dataset could not be built (no shared variants, PLINK failure)."""


def _die(msg: str):
    raise JointDatasetError(msg)


def build_joint_dataset(plink, ref_bfile, user_bfile, needed_ids, work_dir, threads="4", memory="8000",
                        site_policy: str = "wgs", chroms=None, prefix: str = "_joint", cleanup: bool = True):
    """See module docstring. `chroms` restricts both sides (list of '1'..'22');
    None = no restriction. Returns (joint_prefix: Path, user_id_map: dict, info: dict)."""
    if site_policy not in SITE_POLICIES:
        raise ValueError(f"site_policy must be one of {SITE_POLICIES}, got {site_policy!r}")
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    needed_ids = set(needed_ids)

    user_bim_ids = read_bim_ids(str(user_bfile) + ".bim")
    n_user_present = len(needed_ids & user_bim_ids)
    if site_policy == "array":
        needed_ids = needed_ids & user_bim_ids
        if not needed_ids:
            _die("array site policy: the user's genotyped sites share no variant with the requested set")

    extract_list = work_dir / f"{prefix}_needed.txt"
    extract_list.write_text("\n".join(sorted(needed_ids)) + "\n")
    common = ["--allow-extra-chr", "--threads", threads, "--memory", memory]
    chr_args = ["--chr", ",".join(chroms)] if chroms else []

    ref_subset = work_dir / f"{prefix}_ref"
    user_subset = work_dir / f"{prefix}_user"

    def extract(src, out, exclude=None):
        args = ["--bfile", src, "--extract", extract_list, "--keep-allele-order",
                "--make-bed", "--out", out] + chr_args + common
        if exclude is not None:
            args += ["--exclude", exclude]
        return run_plink(plink, args)

    r = extract(ref_bfile, ref_subset)
    if r.returncode != 0 or not Path(str(ref_subset) + ".bed").exists():
        _die(f"could not extract the requested variants from the reference:\n{(r.stderr or r.stdout)[-400:]}")
    r = extract(user_bfile, user_subset)
    if r.returncode != 0 or not Path(str(user_subset) + ".bed").exists():
        _die("the user dataset shares no requested variant with the reference"
             f"{' on chromosomes ' + ','.join(chroms) if chroms else ''}:\n{(r.stderr or r.stdout)[-400:]}")
    user_id_map = rename_user_samples(Path(str(user_subset) + ".fam"))

    merged = work_dir / f"{prefix}_merged"
    merge_args = ["--bfile", ref_subset, "--bmerge", user_subset, "--keep-allele-order",
                  "--make-bed", "--out", merged] + common
    r = run_plink(plink, merge_args)
    missnp = Path(str(merged) + "-merge.missnp")
    n_conflicts = 0
    if r.returncode != 0 and missnp.exists():
        n_conflicts = sum(1 for _ in open(missnp))
        for src, out in ((ref_bfile, ref_subset), (user_bfile, user_subset)):
            r2 = extract(src, out, exclude=missnp)
            if r2.returncode != 0:
                _die(f"PLINK failed re-extracting without .missnp variants:\n{(r2.stderr or r2.stdout)[-400:]}")
        user_id_map = rename_user_samples(Path(str(user_subset) + ".fam"))
        r = run_plink(plink, merge_args)
    if r.returncode != 0 or not Path(str(merged) + ".bed").exists():
        _die(f"PLINK failed merging the user into the reference:\n{(r.stderr or r.stdout)[-600:]}")

    joint = work_dir / prefix
    if site_policy == "wgs":
        # PLINK refuses --fill-missing-a2 alongside --bmerge: second pass.
        r = run_plink(plink, ["--bfile", merged, "--fill-missing-a2", "--keep-allele-order",
                              "--make-bed", "--out", joint] + common)
        if r.returncode != 0 or not Path(str(joint) + ".bed").exists():
            _die(f"PLINK failed filling missing calls as homozygous reference:\n{(r.stderr or r.stdout)[-600:]}")
    else:
        for ext in (".bed", ".bim", ".fam"):
            Path(str(merged) + ext).replace(Path(str(joint) + ext))

    n_variants = sum(1 for _ in open(str(joint) + ".bim"))
    info = {"site_policy": site_policy, "n_requested": len(needed_ids), "n_joint_variants": n_variants,
            "n_user_genotyped": n_user_present, "n_merge_conflicts_excluded": n_conflicts,
            "chromosomes": chroms, "n_samples": len(read_fam_ids(str(joint) + ".fam")),
            "user_bim_ids": user_bim_ids}
    if cleanup:
        for p in (user_subset, merged):
            for ext in (".bed", ".bim", ".fam", ".log", ".nosex", "-merge.missnp"):
                Path(str(p) + ext).unlink(missing_ok=True)
    return joint, user_id_map, info


def remove_joint_dataset(joint_prefix, extra_prefixes=()):
    for p in (joint_prefix, *extra_prefixes):
        for ext in (".bed", ".bim", ".fam", ".log", ".nosex", "-merge.missnp", "_needed.txt"):
            Path(str(p) + ext).unlink(missing_ok=True)


__all__ = ["JointDatasetError", "build_joint_dataset", "remove_joint_dataset", "read_bim_ids", "read_fam_ids",
           "user_covered_chromosomes", "rename_user_samples", "run_plink", "SITE_POLICIES", "AUTOSOMES"]
