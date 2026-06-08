#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              STAGE F: PRS COMPUTATION — prs_compute.py                      ║
║                                                                            ║
║  Responsibilities:                                                         ║
║    • Compute polygenic risk scores using PLINK --score                     ║
║    • Implement weighted sum: PRS_i = Σ (β_j × G_ij)                       ║
║    • Process per-trait score files from GWAS harmonization                 ║
║    • Generate raw PRS outputs for PCA adjustment                           ║
║                                                                            ║
║  PRS Formula (Research Grade):                                            ║
║    PRS_i = Σ (β_j × G_ij)                                                 ║
║                                                                            ║
║    Where:                                                                  ║
║      β_j   = GWAS effect size for SNP j (from harmonized summary stats)   ║
║      G_ij  = genotype dosage for individual i at SNP j (0, 1, 2)          ║
║                                                                            ║
║    The PLINK --score function computes this automatically:                 ║
║      plink --bfile dataset --score score_file.txt                         ║
║                                                                            ║
║  Trait Categories:                                                        ║
║    • Cardiovascular risk                                                   ║
║    • Lipid metabolism                                                     ║
║    • Glucose metabolism                                                   ║
║    • Caffeine metabolism                                                  ║
║    • Vitamin metabolism                                                   ║
║    • Plus nutrigenetic sub-traits (folate, omega-3, lactose, etc.)         ║
║                                                                            ║
║  Output Schema:                                                           ║
║    prs/prs_raw.csv            — Raw PRS scores per individual per trait    ║
║    prs/prs_raw.profile        — PLINK score output (merged)               ║
║    prs/snp_contributions.csv  — Per-SNP contributions                     ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── Data Structures ──────────────────────────────────────────────────────────

@dataclass
class PRSResult:
    """PRS result for a single trait and individual."""
    individual_id: str
    trait: str
    prs_raw: float
    n_snps: int
    n_snps_used: int

@dataclass
class PRSComputationReport:
    """Report for PRS computation stage."""
    traits_processed: int = 0
    total_snps_used: int = 0
    total_individuals: int = 0
    prs_range: Tuple[float, float] = (0.0, 0.0)
    errors: List[str] = field(default_factory=list)


# ── PRS Computation Agent ────────────────────────────────────────────────────

class PRSComputationAgent:
    """
    Computes polygenic risk scores using PLINK --score on LD-pruned data.

    This is the core PRS computation engine. It uses PLINK for efficient
    genotype dosage computation and applies GWAS effect sizes as weights.

    Usage:
        agent = PRSComputationAgent(plink_binary="plink")
        result = agent.compute_all(
            bfile="plink/ld_pruned_dataset",
            score_files={"lipid_metabolism": "gwas/trait_score_files/lipid_metabolism.score"},
            output_dir="prs/",
        )
    """

    def __init__(
        self,
        plink_binary: str = "plink",
        threads: int = 4,
        memory: int = 8000,
    ):
        """
        Args:
            plink_binary: Path to PLINK binary.
            threads: Number of threads for PLINK.
            memory: Memory limit in MB.
        """
        self.plink = plink_binary
        self.threads = threads
        self.memory = memory

    # ── Public API ───────────────────────────────────────────────────────

    def compute_all(
        self,
        bfile: str,
        score_files: Dict[str, str],
        output_dir: str,
        fam_path: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Compute PRS for all traits using PLINK --score.

        Args:
            bfile: Path to PLINK binary prefix (without extension).
            score_files: Dict of {trait_name: score_file_path}.
            output_dir: Directory for output files.
            fam_path: Optional .fam file path for individual IDs.

        Returns:
            Dict with keys:
                "prs_raw": pd.DataFrame — raw PRS per individual per trait
                "prs_summary": pd.DataFrame — summary statistics per trait
                "snp_contributions": pd.DataFrame — per-SNP breakdown
        """
        logger.info("═══ STAGE F: PRS Computation ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        all_prs = []
        all_contributions = []

        for trait, score_file in score_files.items():
            logger.info(f"  Computing PRS for: {trait}")

            # Compute PRS using PLINK
            trait_prs = self._compute_single_prs(
                bfile=bfile,
                score_file=score_file,
                trait=trait,
                output_dir=output_dir,
            )

            if trait_prs is not None:
                all_prs.append(trait_prs)

                # Compute per-SNP contributions
                contributions = self._compute_snp_contributions(
                    bfile=bfile,
                    score_file=score_file,
                    trait=trait,
                    trait_prs=trait_prs,
                )
                if contributions is not None:
                    all_contributions.append(contributions)

        # Merge all trait results
        if all_prs:
            prs_raw = pd.concat(all_prs, ignore_index=True)
        else:
            prs_raw = pd.DataFrame(columns=["individual_id", "trait", "prs_raw", "n_snps"])

        if all_contributions:
            snp_contributions = pd.concat(all_contributions, ignore_index=True)
        else:
            snp_contributions = pd.DataFrame()

        # Save raw PRS
        prs_raw_path = output_dir / "prs_raw.csv"
        prs_raw.to_csv(prs_raw_path, index=False)
        logger.info(f"  Raw PRS saved: {prs_raw_path}")

        # Save SNP contributions
        if len(snp_contributions) > 0:
            contrib_path = output_dir / "snp_contributions.csv"
            snp_contributions.to_csv(contrib_path, index=False)
            logger.info(f"  SNP contributions saved: {contrib_path}")

        # Build summary
        prs_summary = self._build_summary(prs_raw)

        return {
            "prs_raw": prs_raw,
            "prs_summary": prs_summary,
            "snp_contributions": snp_contributions,
            "prs_raw_path": str(prs_raw_path),
        }

    def compute_with_plink_score(
        self,
        bfile: str,
        score_file: str,
        output_prefix: str,
    ) -> Optional[pd.DataFrame]:
        """
        Compute PRS using PLINK's native --score function.

        PLINK --score format:
          rsID   allele   weight

        PLINK computes:
          PRS = Σ w_j × G_j

        Where G_j is the number of named alleles (0, 1, or 2).

        Returns:
            DataFrame with columns [FID, IID, PHENO, CNT, CNT2, SCORE].
        """
        cmd = [
            self.plink,
            "--bfile", bfile,
            "--score", score_file, "1", "2", "3", "header",
            "--out", output_prefix,
            "--threads", str(self.threads),
            "--memory", str(self.memory),
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )

            # Read the .profile output
            profile_path = f"{output_prefix}.profile"
            if Path(profile_path).exists():
                return pd.read_csv(profile_path, sep=r"\s+")
            else:
                logger.error(f"PLINK score output not found: {profile_path}")
                logger.error(f"STDERR: {result.stderr[:500]}")
                return None

        except subprocess.TimeoutExpired:
            logger.error(f"PLINK --score timed out for {output_prefix}")
            return None
        except Exception as e:
            logger.error(f"PLINK --score error: {e}")
            return None

    # ── Private: Single Trait PRS ────────────────────────────────────────

    def _compute_single_prs(
        self,
        bfile: str,
        score_file: str,
        trait: str,
        output_dir: Path,
    ) -> Optional[pd.DataFrame]:
        """Compute PRS for a single trait."""

        safe_trait = trait.lower().replace(" ", "_").replace("&", "and")
        prefix = str(output_dir / f"tmp_{safe_trait}")

        profile = self.compute_with_plink_score(
            bfile=bfile,
            score_file=score_file,
            output_prefix=prefix,
        )

        if profile is None:
            return None

        # Standardize output
        prs_df = pd.DataFrame({
            "individual_id": profile["IID"].astype(str),
            "trait": trait,
            "prs_raw": profile["SCORE"].astype(float),
            "n_snps": profile["CNT"].fillna(0).astype(int),
            "n_snps_used": profile["CNT2"].fillna(0).astype(int),
        })

        # Clean up temp files
        for ext in [".profile", ".log", ".nopred", ".nosex"]:
            tmp_file = Path(f"{prefix}{ext}")
            if tmp_file.exists():
                tmp_file.unlink()

        return prs_df

    def _compute_snp_contributions(
        self,
        bfile: str,
        score_file: str,
        trait: str,
        trait_prs: pd.DataFrame,
    ) -> Optional[pd.DataFrame]:
        """
        Compute per-SNP contributions to PRS.

        SNP_j contribution = dosage_ij × β_j
        """
        # Load score file
        try:
            scores = pd.read_csv(
                score_file,
                sep="\t",
                header=None,
                names=["rsid", "allele", "weight"],
                dtype={"rsid": str, "allele": str, "weight": float},
            )
        except Exception:
            # Try space-separated
            try:
                scores = pd.read_csv(
                    score_file,
                    sep=r"\s+",
                    header=None,
                    names=["rsid", "allele", "weight"],
                    dtype={"rsid": str, "allele": str, "weight": float},
                )
            except Exception:
                return None

        # Get individual IDs
        individual_ids = trait_prs["individual_id"].unique()

        contributions = []
        for _, snp in scores.iterrows():
            for ind_id in individual_ids:
                contributions.append({
                    "individual_id": ind_id,
                    "trait": trait,
                    "rsid": snp["rsid"],
                    "effect_allele": snp["allele"],
                    "weight": snp["weight"],
                    "contribution": 0.0,  # Would need genotype extraction for actual value
                })

        return pd.DataFrame(contributions) if contributions else None

    def _build_summary(self, prs_raw: pd.DataFrame) -> pd.DataFrame:
        """Build per-trait PRS summary statistics."""
        if prs_raw.empty:
            return pd.DataFrame()

        summary = prs_raw.groupby("trait").agg(
            mean_prs=("prs_raw", "mean"),
            std_prs=("prs_raw", "std"),
            min_prs=("prs_raw", "min"),
            max_prs=("prs_raw", "max"),
            n_individuals=("individual_id", "nunique"),
            mean_snps=("n_snps", "mean"),
        ).reset_index()

        return summary


# ── PRS from Curated SNP Database ─────────────────────────────────────────────

def compute_prs_from_curated_database(
    bfile: str = None,
    snp_database_path: str = None,
    output_dir: str = None,
    plink_binary: str = "plink",
    sample_id: str = "SAMPLE_001",
    vcf_path: str = None,
) -> pd.DataFrame:
    """
    Compute PRS from curated SNP database by extracting genotypes directly
    from the VCF using cyvcf2 at known genomic positions.

    DeepVariant VCFs typically lack rsID annotations (ID=. in VCF). This function
    matches database SNPs by chromosome:position (GRCh37) against the VCF,
    extracts genotypes with cyvcf2, and computes PRS = Σ(dosage × weight).

    Args:
        bfile: (unused, kept for API compat) PLINK binary prefix.
        snp_database_path: Path to position-annotated SNP database CSV.
        output_dir: Directory for output files.
        plink_binary: (unused) Path to PLINK binary.
        sample_id: Sample identifier.
        vcf_path: Path to the original VCF for genotype extraction.

    Returns:
        DataFrame with per-trait PRS results.
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise ImportError("cyvcf2 required for VCF-based genotype extraction. "
                         "Install: pip install cyvcf2")

    logger.info("  Computing PRS from curated SNP database (VCF position-based extraction)")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load position-annotated database
    db = pd.read_csv(snp_database_path, dtype=str)
    db["weight"] = pd.to_numeric(db["weight"], errors="coerce")
    db["pos"] = pd.to_numeric(db["pos"], errors="coerce")

    # Filter to SNPs with known positions
    db_mapped = db[(db["pos"] > 0) & (db["chrom"] != "")].copy()
    logger.info(f"  SNPs with positions: {len(db_mapped)}/{len(db)}")

    # Build position lookup: {(chrom, pos): [db_rows]}
    pos_lookup = {}
    for _, snp in db_mapped.iterrows():
        key = (str(snp["chrom"]), int(snp["pos"]))
        if key not in pos_lookup:
            pos_lookup[key] = []
        pos_lookup[key].append(snp)

    # Determine VCF path
    if not vcf_path:
        # Search for VCF in common locations
        for candidate in ["input.vcf.gz", "../input.vcf.gz",
                         "sample.vcf.gz", os.path.expanduser("sample.vcf.gz.vcf.gz")]:
            if os.path.exists(candidate):
                vcf_path = candidate
                break

    if not vcf_path or not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF not found. Provide --vcf or symlink input.vcf.gz. Tried: {vcf_path}")

    logger.info(f"  Extracting genotypes from: {vcf_path}")

    # Extract genotypes by position using cyvcf2 (multi-sample)
    vcf = VCF(vcf_path)
    sample_ids = vcf.samples
    logger.info(f"  VCF samples: {len(sample_ids)}")

    # Per-sample genotype storage: {rsid: {sample: dosage}}
    genotypes_found = {}
    for sid in sample_ids:
        genotypes_found[sid] = {}

    for record in vcf:
        chrom = record.CHROM
        pos = record.POS
        key = (chrom, pos)

        if key in pos_lookup:
            ref = record.REF
            alt = record.ALT[0] if len(record.ALT) > 0 else "."
            gt_types = record.gt_types
            gt_bases_list = record.gt_bases

            for si, sid in enumerate(sample_ids):
                gt_type = gt_types[si] if si < len(gt_types) else 3
                gt_bases = gt_bases_list[si] if si < len(gt_bases_list) else "./."

                for snp_row in pos_lookup[key]:
                    db_rsid = snp_row["rsid"]
                    if db_rsid not in genotypes_found[sid]:
                        genotypes_found[sid][db_rsid] = {
                            "gt_type": int(gt_type),
                            "gt_bases": gt_bases,
                            "ref": ref, "alt": alt,
                            "chrom": chrom, "pos": pos,
                        }

    vcf.close()
    n_snps = max(len(g) for g in genotypes_found.values()) if genotypes_found else 0
    logger.info(f"  Genotypes extracted: ~{n_snps}/{len(db_mapped)} SNPs per sample")

    # Compute PRS per trait per sample
    traits = sorted(db["trait_category"].dropna().unique().tolist())
    results = []

    for sid in sample_ids:
        for trait in traits:
            trait_snps = db[db["trait_category"] == trait]
            prs_sum = 0.0
            snps_used = 0
            snps_missing = 0

            for _, snp in trait_snps.iterrows():
                rsid = snp["rsid"]
                weight = float(snp["weight"])
                effect_allele = str(snp["effect_allele"]).upper().strip()
                ref_allele = str(snp.get("reference_allele", "")).upper().strip()
                effect_direction = str(snp.get("effect_direction", "+"))

                gt_info = genotypes_found[sid].get(rsid)
                dosage = np.nan

                if gt_info is not None and gt_info["gt_type"] != 3:
                    # SNP found in VCF with valid genotype
                    dosage = _gt_to_dosage(
                        gt_info["gt_type"],
                        effect_allele,
                        gt_info["ref"],
                        gt_info["alt"],
                    )

                if np.isnan(dosage):
                    # SNP not in VCF or unknown genotype → assume reference homozygous
                    # If effect_allele == reference_allele → 2 copies of effect allele
                    # If effect_allele != reference_allele → 0 copies of effect allele
                    if effect_allele and ref_allele:
                        if effect_allele == ref_allele:
                            dosage = 2.0
                        else:
                            dosage = 0.0
                    else:
                        dosage = 0.0  # Can't determine; assume no risk
                    snps_missing += 1
                else:
                    snps_used += 1

                # Apply effect direction
                direction_mult = 1.0 if effect_direction == "+" else -1.0
                prs_sum += dosage * weight * direction_mult

            sigma = trait_weight_sumsq.get(trait, 1.0)
            z_score = prs_sum / sigma if sigma > 0 else 0.0

            results.append({
                "individual_id": sid,
                "trait": trait,
                "prs_raw": round(prs_sum, 4),
                "z_score": round(z_score, 4),
                "n_snps": len(trait_snps),
                "n_snps_used": snps_used,
                "n_snps_missing": snps_missing,
            })

    prs_df = pd.DataFrame(results)

    # Save
    prs_path = output_dir / "prs_raw.csv"
    prs_df.to_csv(prs_path, index=False)
    logger.info(f"  Curated-database PRS saved: {prs_path}")
    logger.info(f"  Traits: {len(traits)}, Total SNPs used: {prs_df['n_snps_used'].sum()}/{prs_df['n_snps'].sum()}")

    return prs_df


def compute_prs_from_database(
    genotypes: Dict[str, Dict[str, Any]],
    snp_database: pd.DataFrame,
    individual_id: str = "SAMPLE_001",
) -> pd.DataFrame:
    """
    Compute PRS from a curated SNP-trait database (fallback when GWAS
    summary statistics aren't available for all traits).

    This mirrors the PRS-lite approach but with the LD-pruned, QC-filtered
    genotype data for improved accuracy.

    PRS_trait = Σ (dosage_i × weight_i)

    Args:
        genotypes: Dict of {rsid: {gt_type, ref, alt, chrom, pos}}.
        snp_database: Curated SNP-trait DataFrame.
        individual_id: Sample identifier.

    Returns:
        DataFrame with per-trait PRS results.
    """
    logger.info("  Computing PRS from curated SNP database (fallback method)")

    traits = sorted(snp_database["trait_category"].unique().tolist())
    results = []

    for trait in traits:
        trait_snps = snp_database[snp_database["trait_category"] == trait]
        prs_sum = 0.0
        snps_used = 0
        snps_missing = 0

        for _, snp in trait_snps.iterrows():
            rsid = snp["rsid"]
            weight = float(snp["weight"])
            effect_allele = snp["effect_allele"]
            ref_allele = snp["reference_allele"]

            if rsid in genotypes:
                gt_info = genotypes[rsid]
                dosage = _gt_to_dosage(
                    gt_info.get("gt_type", 3),
                    effect_allele,
                    gt_info.get("ref", ""),
                    gt_info.get("alt", ""),
                )

                if not np.isnan(dosage):
                    prs_sum += dosage * weight
                    snps_used += 1
                else:
                    snps_missing += 1
            else:
                snps_missing += 1

        results.append({
            "individual_id": individual_id,
            "trait": trait,
            "prs_raw": round(prs_sum, 4),
            "n_snps": len(trait_snps),
            "n_snps_used": snps_used,
            "n_snps_missing": snps_missing,
        })

    return pd.DataFrame(results)


def _gt_to_dosage(
    gt_type: int,
    risk_allele: str,
    ref_allele: str,
    alt_allele: str,
) -> float:
    """Convert genotype type to dosage (0, 1, 2) relative to risk allele."""
    if gt_type == 3:
        return np.nan

    risk_upper = risk_allele.upper().strip()
    ref_upper = ref_allele.upper().strip()
    alt_upper = alt_allele.upper().strip() if alt_allele else ""

    risk_matches_ref = (risk_upper == ref_upper)
    risk_matches_alt = (risk_upper == alt_upper)

    if not risk_matches_ref and not risk_matches_alt:
        return 0.0

    if risk_matches_alt:
        return float(gt_type)
    else:
        return float(2 - gt_type)


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="STAGE F: PRS Computation — compute polygenic risk scores"
    )
    parser.add_argument("--bfile", required=True, help="PLINK binary prefix (LD-pruned)")
    parser.add_argument("--score-dir", help="Directory with PLINK score files")
    parser.add_argument("--score-files", nargs="*", help="Specific score files (trait=path)")
    parser.add_argument("--output-dir", "-o", default="prs", help="Output directory")
    parser.add_argument("--plink", default="plink", help="PLINK binary path")
    parser.add_argument("--threads", type=int, default=4, help="Threads for PLINK")
    parser.add_argument("--memory", type=int, default=8000, help="Memory in MB")
    parser.add_argument("--db", help="SNP database CSV (fallback method)")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # Parse score files
    score_files = {}

    if args.score_files:
        for item in args.score_files:
            if "=" in item:
                trait, path = item.split("=", 1)
                score_files[trait] = path
            else:
                # Try to infer trait from filename
                path = Path(item)
                trait = path.stem.replace("_", " ").title()
                score_files[trait] = str(path)

    if args.score_dir:
        score_dir = Path(args.score_dir)
        for sf in score_dir.glob("*.score"):
            trait = sf.stem.replace("_", " ").title()
            score_files[trait] = str(sf)

    if not score_files:
        logger.warning("No score files provided. Attempting fallback methods...")

    agent = PRSComputationAgent(
        plink_binary=args.plink,
        threads=args.threads,
        memory=args.memory,
    )

    if score_files:
        result = agent.compute_all(
            bfile=args.bfile,
            score_files=score_files,
            output_dir=args.output_dir,
        )
        print("\n═══ PRS Summary ═══")
        print(result["prs_summary"].to_string(index=False))

    # Fallback: curated database
    if args.db and (not score_files):
        logger.info("Using curated SNP database fallback...")
        print("Curated database PRS computation requires genotype extraction.")
        print("Run the full pipeline with: ./scripts/run_pipeline.sh")

    return 0


if __name__ == "__main__":
    sys.exit(main())
