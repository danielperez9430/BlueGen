#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   UNCERTAINTY PROPAGATION ENGINE — scripts/14_uncertainty_propagation.py     ║
║                                                                            ║
║   Propagates uncertainty through the entire PRS pipeline:                  ║
║                                                                            ║
║   Layer 1: Genotype uncertainty (GQ/DP from VCF)                           ║
║   Layer 2: Ancestry uncertainty (PCA softmax probabilities)                ║
║   Layer 3: GWAS effect uncertainty (SE from summary stats or evidence)     ║
║                                                                            ║
║   Propagation rule:                                                        ║
║     Var(PRS_i) = Σ_j [β_j² × Var(G_ij)] + Σ_j [G_ij² × Var(β_j)]         ║
║                                                                            ║
║   Total uncertainty = sqrt(Var_geno + Var_ancestry + Var_effect)           ║
║                                                                            ║
║   Every output value includes:                                             ║
║     • point estimate                                                       ║
║     • 95% confidence interval                                              ║
║     • per-layer uncertainty decomposition                                  ║
║     • uncertainty_score (normalized 0–1)                                   ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# Data Structures
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class UncertaintyDecomposition:
    """Breakdown of uncertainty sources."""
    total_variance: float
    genotype_variance: float
    ancestry_variance: float
    effect_variance: float
    genotype_fraction: float       # % of total variance from genotypes
    ancestry_fraction: float
    effect_fraction: float


@dataclass
class PRSWithUncertainty:
    """PRS value with full uncertainty quantification."""
    individual_id: str
    trait: str
    prs_point_estimate: float
    prs_std_error: float
    confidence_interval_95: Tuple[float, float]
    confidence_interval_68: Tuple[float, float]
    uncertainty_score: float              # 0–1, normalized
    decomposition: UncertaintyDecomposition
    n_snps_with_genotype: int
    n_snps_with_effect_se: int


@dataclass
class UncertaintyReport:
    """Complete uncertainty propagation report."""
    results: List[PRSWithUncertainty]
    global_uncertainty_score: float
    method: str
    genotype_quality_summary: Dict[str, float]
    ancestry_entropy: float
    gwas_evidence_summary: Dict[str, Any]


# ═══════════════════════════════════════════════════════════════════════════════
# Uncertainty Propagation Engine
# ═══════════════════════════════════════════════════════════════════════════════

class UncertaintyPropagationEngine:
    """
    Propagates uncertainty through three layers into final PRS outputs.

    Layer 1 — Genotype uncertainty:
      Var(G) derived from Phred-scaled genotype likelihoods (PL field in VCF)
      or GQ (Genotype Quality). For missing genotypes, use max uncertainty.

    Layer 2 — Ancestry uncertainty:
      Entropy of ancestry probability distribution → uncertainty in
      population-specific PRS calibration parameters (mu, sigma).

    Layer 3 — GWAS effect uncertainty:
      From SE in GWAS summary statistics, or estimated from evidence level
      (A: SE ≈ |β|/5, B: SE ≈ |β|/3, C: SE ≈ |β|/2, D: SE ≈ |β|).

    Propagation:
      Var(PRS) = Σ β²·Var(G) + Σ G²·Var(β)
      Total SE = sqrt(Var_geno + Var_ancestry + Var_effect)

    Usage:
        engine = UncertaintyPropagationEngine(vcf_path="input.vcf.gz")
        result = engine.propagate(
            prs_data=prs_df,
            ancestry_path="pca/ancestry_inference.json",
            snp_database_path="data/snp_database_annotated.csv",
        )
    """

    # Evidence level → approximate SE/|β| ratio
    # Calibrated for modern GWAS precision (v1.2.0 — reduced conservatism)
    EVIDENCE_SE_RATIO = {
        "A": 0.10,    # GWAS p < 5e-8, n>50k: tight CI (SE ≈ |β|/10)
        "B": 0.20,    # Replicated: good precision (SE ≈ |β|/5)
        "C": 0.35,    # Single study: moderate (SE ≈ |β|/3)
        "D": 0.55,    # Mechanistic: wide but bounded (SE ≈ |β|/2)
    }

    # GQ → genotype uncertainty mapping
    GQ_UNCERTAINTY = {
        # GQ range → P(genotype error)
        (0, 10):   0.50,
        (10, 20):  0.25,
        (20, 30):  0.10,
        (30, 40):  0.05,
        (40, 100): 0.01,
        (100, 999): 0.005,
    }

    def __init__(self, vcf_path: Optional[str] = None):
        self.vcf_path = vcf_path
        self._genotype_uncertainties: Dict[str, float] = {}

    # ── Public API ───────────────────────────────────────────────────────

    def propagate(
        self,
        prs_data: pd.DataFrame,
        ancestry_path: str,
        snp_database_path: str,
        output_dir: str,
        consistency_check_path: Optional[str] = None,
        reference_distributions_path: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Propagate uncertainty through all layers.

        Args:
            prs_data: PRS DataFrame with [individual_id, trait, prs_raw, prs_adjusted].
            ancestry_path: Ancestry inference JSON.
            snp_database_path: Position-annotated SNP database CSV.
            output_dir: Output directory.
            consistency_check_path: GWAS consistency check JSON (optional).
            reference_distributions_path: Population reference distributions JSON (for population SD).

        Returns:
            Dict with uncertainty-quantified PRS and report.
        """
        logger.info("═══ Uncertainty Propagation Engine ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # ── Load population SD per trait for uncertainty denominator ────────
        trait_population_sd = {}
        if reference_distributions_path is None:
            reference_distributions_path = "reference/population_distributions/reference_distributions.json"
        if os.path.exists(reference_distributions_path):
            try:
                with open(reference_distributions_path) as fh:
                    ref_dists = json.load(fh)
                for trait, pops in ref_dists.get("distributions", {}).items():
                    eur = pops.get("EUR", {})
                    sd = eur.get("std", 0)
                    if sd and sd > 0:
                        trait_population_sd[trait.lower()] = sd
                logger.info(f"  Population SD loaded for {len(trait_population_sd)} traits")
            except Exception as e:
                logger.warning(f"  Could not load population distributions: {e}")
        else:
            logger.info(f"  No population distributions at {reference_distributions_path}")

        # ── Layer 1: Genotype uncertainty ────────────────────────────────
        geno_uncertainty = self._compute_genotype_uncertainty(snp_database_path)
        geno_summary = self._summarize_genotype_quality(geno_uncertainty)
        logger.info(f"  Genotype uncertainty: mean={geno_summary.get('mean_uncertainty', 0):.4f}")

        # ── Layer 2: Ancestry uncertainty ─────────────────────────────────
        with open(ancestry_path) as fh:
            ancestry = json.load(fh)
        ancestry_probs = ancestry.get("summary", {}).get("all_probabilities", {})
        ancestry_entropy = self._compute_ancestry_entropy(ancestry_probs)
        ancestry_uncertainty = self._compute_ancestry_uncertainty(ancestry_probs)
        logger.info(f"  Ancestry entropy: {ancestry_entropy:.4f}")
        logger.info(f"  Ancestry uncertainty (sigma_mu): {ancestry_uncertainty:.4f}")

        # ── Layer 3: GWAS effect uncertainty ─────────────────────────────
        effect_uncertainty = self._compute_effect_uncertainty(snp_database_path)
        logger.info(f"  Effect uncertainty: {len(effect_uncertainty)} SNPs quantified")

        # ── Load SNP database for weights ─────────────────────────────────
        db = pd.read_csv(snp_database_path, dtype=str)
        db["weight"] = pd.to_numeric(db["weight"], errors="coerce")

        # ── Consistency check downgrade ───────────────────────────────────
        consistency_downgrade = 0.0
        if consistency_check_path and os.path.exists(consistency_check_path):
            with open(consistency_check_path) as fh:
                cc = json.load(fh)
            consistency_downgrade = cc.get("confidence_downgrade", 0.0)
            logger.info(f"  Consistency downgrade: {consistency_downgrade:.0%}")

        # ── Propagate per trait ──────────────────────────────────────────
        results = []
        prs_col = "prs_adjusted" if "prs_adjusted" in prs_data.columns else "prs_raw"

        for _, row in prs_data.iterrows():
            trait = row["trait"]
            iid = row.get("individual_id", "SAMPLE_001")
            prs_val = float(row[prs_col])

            # Get trait-specific SNPs
            trait_snps = db[db["trait_category"] == trait]

            # Compute uncertainty components
            geno_var = self._propagate_genotype_variance(trait_snps, geno_uncertainty)
            effect_var = self._propagate_effect_variance(trait_snps, effect_uncertainty)
            ancestry_var = ancestry_uncertainty ** 2  # Affects mu/sigma calibration

            total_var = geno_var + effect_var + ancestry_var
            total_se = np.sqrt(max(total_var, 1e-10))

            # Apply consistency downgrade
            total_se *= (1.0 + consistency_downgrade)
            total_var = total_se ** 2

            # Confidence intervals
            ci_95 = (prs_val - 1.96 * total_se, prs_val + 1.96 * total_se)
            ci_68 = (prs_val - total_se, prs_val + total_se)

            # Uncertainty score (0–1)
            # Denominator: population-based (SD × 2 captures ~95% of population range)
            # Much more stable than the old abs(prs_val) * 0.5 which saturated for low-PRS individuals
            population_sd = trait_population_sd.get(trait.lower(), 0.25)
            max_plausible_se = max(population_sd * 2.0, 0.1)
            uncertainty_score = min(total_se / max(max_plausible_se, 1e-6), 1.0)

            # Decomposition
            decomp = UncertaintyDecomposition(
                total_variance=round(total_var, 6),
                genotype_variance=round(geno_var, 6),
                ancestry_variance=round(ancestry_var, 6),
                effect_variance=round(effect_var, 6),
                genotype_fraction=round(geno_var / max(total_var, 1e-10), 4),
                ancestry_fraction=round(ancestry_var / max(total_var, 1e-10), 4),
                effect_fraction=round(effect_var / max(total_var, 1e-10), 4),
            )

            n_geno = sum(1 for _, s in trait_snps.iterrows()
                        if s["rsid"] in geno_uncertainty)
            n_se = sum(1 for _, s in trait_snps.iterrows()
                      if s["rsid"] in effect_uncertainty)

            results.append(PRSWithUncertainty(
                individual_id=iid,
                trait=trait,
                prs_point_estimate=round(prs_val, 4),
                prs_std_error=round(total_se, 4),
                confidence_interval_95=(round(ci_95[0], 4), round(ci_95[1], 4)),
                confidence_interval_68=(round(ci_68[0], 4), round(ci_68[1], 4)),
                uncertainty_score=round(uncertainty_score, 4),
                decomposition=decomp,
                n_snps_with_genotype=n_geno,
                n_snps_with_effect_se=n_se,
            ))

        # Global uncertainty score
        global_score = np.mean([r.uncertainty_score for r in results]) if results else 0.0

        # Build report
        report = UncertaintyReport(
            results=results,
            global_uncertainty_score=round(global_score, 4),
            method="three_layer_variance_propagation",
            genotype_quality_summary=geno_summary,
            ancestry_entropy=round(ancestry_entropy, 4),
            gwas_evidence_summary=self._summarize_evidence(effect_uncertainty),
        )

        # ── Save outputs ──────────────────────────────────────────────────
        # Uncertainty-augmented PRS table
        output_rows = []
        for r in results:
            output_rows.append({
                "individual_id": r.individual_id,
                "trait": r.trait,
                "prs": r.prs_point_estimate,
                "prs_se": r.prs_std_error,
                "ci_95_lower": r.confidence_interval_95[0],
                "ci_95_upper": r.confidence_interval_95[1],
                "uncertainty_score": r.uncertainty_score,
                "genotype_variance": r.decomposition.genotype_variance,
                "ancestry_variance": r.decomposition.ancestry_variance,
                "effect_variance": r.decomposition.effect_variance,
                "genotype_fraction": r.decomposition.genotype_fraction,
            })

        uncertainty_df = pd.DataFrame(output_rows)
        uncer_path = output_dir / "prs_uncertainty.csv"
        uncertainty_df.to_csv(uncer_path, index=False)
        logger.info(f"  Uncertainty PRS: {uncer_path}")

        # JSON report
        report_path = output_dir / "uncertainty_report.json"
        with open(report_path, "w") as fh:
            json.dump({
                "global_uncertainty_score": report.global_uncertainty_score,
                "method": report.method,
                "genotype_quality_summary": report.genotype_quality_summary,
                "ancestry_entropy": report.ancestry_entropy,
                "results": [asdict(r) for r in results],
            }, fh, indent=2, default=str)
        logger.info(f"  Uncertainty report: {report_path}")

        return {
            "prs_uncertainty": uncertainty_df,
            "report": report,
            "uncertainty_path": str(uncer_path),
            "report_path": str(report_path),
        }

    # ── Private: Genotype Uncertainty ────────────────────────────────────

    def _compute_genotype_uncertainty(
        self, snp_database_path: str
    ) -> Dict[str, float]:
        """
        Extract genotype uncertainty from VCF GQ field.
        Returns {rsid: Var(G)} for each SNP.
        """
        uncertainties = {}

        if not self.vcf_path or not os.path.exists(self.vcf_path):
            logger.warning("  No VCF available — using default genotype uncertainty")
            return uncertainties

        try:
            from cyvcf2 import VCF
            db = pd.read_csv(snp_database_path, dtype=str)
            db["pos"] = pd.to_numeric(db["pos"], errors="coerce")

            pos_set = set()
            for _, s in db.iterrows():
                if s["pos"] > 0:
                    pos_set.add((str(s["chrom"]), int(s["pos"]), s["rsid"]))

            vcf = VCF(self.vcf_path)
            for record in vcf:
                key = (record.CHROM, record.POS)
                for chrom, pos, rsid in list(pos_set):
                    if key == (chrom, pos):
                        gq = record.gt_quals[0] if len(record.gt_quals) > 0 else 0
                        gt_type = record.gt_types[0] if len(record.gt_types) > 0 else 3

                        # Convert GQ to genotype error probability
                        p_error = self._gq_to_error_prob(gq)

                        # Var(G) for dosage (0/1/2) given error probability
                        # Under Hardy-Weinberg, but simplified:
                        # Var(G|observed) ≈ 2 * p_error * (1 - p_error)
                        # For unknown genotypes, use max uncertainty
                        if gt_type == 3:
                            var_g = 0.5  # Maximum uncertainty
                        else:
                            var_g = 2.0 * p_error * (1.0 - p_error)

                        uncertainties[rsid] = var_g

            vcf.close()

        except ImportError:
            logger.warning("  cyvcf2 not available — genotype uncertainty estimated from defaults")
        except Exception as e:
            logger.warning(f"  VCF parsing error: {e}")

        # For SNPs not in VCF (reference homozygous), uncertainty is low
        for _, s in db.iterrows():
            rsid = s["rsid"]
            if rsid not in uncertainties:
                uncertainties[rsid] = 0.001  # Very low uncertainty for ref calls

        return uncertainties

    def _gq_to_error_prob(self, gq: float) -> float:
        """Convert Genotype Quality (Phred scale) to error probability."""
        if gq <= 0:
            return 0.50
        p = 10 ** (-gq / 10.0)
        return min(p, 0.50)

    def _summarize_genotype_quality(
        self, uncertainties: Dict[str, float]
    ) -> Dict[str, float]:
        """Summarize genotype quality statistics."""
        if not uncertainties:
            return {"mean_uncertainty": 0.01, "n_snps": 0}

        vals = list(uncertainties.values())
        return {
            "mean_uncertainty": round(np.mean(vals), 6),
            "median_uncertainty": round(np.median(vals), 6),
            "max_uncertainty": round(np.max(vals), 6),
            "min_uncertainty": round(np.min(vals), 6),
            "n_snps_high_quality": sum(1 for v in vals if v < 0.01),
            "n_snps_low_quality": sum(1 for v in vals if v > 0.1),
            "n_snps_total": len(vals),
        }

    # ── Private: Ancestry Uncertainty ────────────────────────────────────

    @staticmethod
    def _compute_ancestry_entropy(probs: Dict[str, float]) -> float:
        """Shannon entropy of ancestry probability distribution."""
        vals = [p for p in probs.values() if p > 0]
        if not vals:
            return 0.0
        return -sum(p * np.log(p) for p in vals)

    @staticmethod
    def _compute_ancestry_uncertainty(probs: Dict[str, float]) -> float:
        """
        Compute ancestry-induced uncertainty in PRS calibration.
        Higher entropy → less certain population parameters → wider CIs.
        """
        entropy = UncertaintyPropagationEngine._compute_ancestry_entropy(probs)
        # Map entropy [0, ~1.6] to uncertainty [0, 0.5]
        max_entropy = np.log(5)  # 5 populations
        normalized = entropy / max_entropy if max_entropy > 0 else 0
        # Scale to PRS units
        return normalized * 0.25

    # ── Private: Effect Uncertainty ──────────────────────────────────────

    def _compute_effect_uncertainty(
        self, snp_database_path: str
    ) -> Dict[str, float]:
        """
        Estimate GWAS effect size uncertainty from evidence levels.
        Returns {rsid: Var(β)}.

        For formal GWAS with SE column, uses actual SE².
        For curated database, estimates from evidence level.
        """
        uncertainties = {}

        try:
            db = pd.read_csv(snp_database_path, dtype=str)
            db["weight"] = pd.to_numeric(db["weight"], errors="coerce")

            for _, snp in db.iterrows():
                rsid = snp["rsid"]
                beta = abs(float(snp["weight"])) if not pd.isna(snp["weight"]) else 0.0

                # Check for actual SE column
                if "se" in snp.index and not pd.isna(snp.get("se")):
                    se = float(snp["se"])
                    uncertainties[rsid] = se ** 2
                else:
                    # Estimate from evidence level
                    evidence = snp.get("evidence_level", "D")
                    se_ratio = self.EVIDENCE_SE_RATIO.get(evidence, 0.50)
                    se = beta * se_ratio
                    uncertainties[rsid] = max(se ** 2, 0.0001)  # Minimum variance

        except Exception as e:
            logger.warning(f"  Effect uncertainty estimation error: {e}")

        return uncertainties

    def _summarize_evidence(
        self, effect_uncertainty: Dict[str, float]
    ) -> Dict[str, Any]:
        """Summarize GWAS evidence quality."""
        if not effect_uncertainty:
            return {}

        se_vals = [np.sqrt(v) for v in effect_uncertainty.values()]
        return {
            "mean_se": round(np.mean(se_vals), 6),
            "median_se": round(np.median(se_vals), 6),
            "n_snps_high_precision": sum(1 for v in se_vals if v < 0.05),
            "n_snps_low_precision": sum(1 for v in se_vals if v > 0.2),
            "n_snps_total": len(se_vals),
        }

    # ── Private: Variance Propagation ────────────────────────────────────

    def _propagate_genotype_variance(
        self,
        trait_snps: pd.DataFrame,
        geno_uncertainty: Dict[str, float],
    ) -> float:
        """Var_geno(PRS) = Σ β² × Var(G)."""
        total = 0.0
        for _, snp in trait_snps.iterrows():
            rsid = snp["rsid"]
            beta = abs(float(snp["weight"])) if not pd.isna(snp.get("weight")) else 0.0
            var_g = geno_uncertainty.get(rsid, 0.01)  # Default: moderate uncertainty
            total += (beta ** 2) * var_g
        return total

    def _propagate_effect_variance(
        self,
        trait_snps: pd.DataFrame,
        effect_uncertainty: Dict[str, float],
    ) -> float:
        """Var_effect(PRS) = Σ G² × Var(β).
        Uses HWE-based E[G²] from MAF instead of hardcoded 1.0."""
        total = 0.0
        for _, snp in trait_snps.iterrows():
            rsid = snp["rsid"]
            var_beta = effect_uncertainty.get(rsid, 0.01)
            # E[G²] = Var(G) + E[G]² = 2·MAF·(1-MAF) + 4·MAF² under HWE
            maf = self._get_maf(snp)
            expected_g2 = 2 * maf * (1 - maf) + 4 * maf**2
            total += expected_g2 * var_beta
        return total

    @staticmethod
    def _get_maf(snp: pd.Series) -> float:
        """Get minor allele frequency for a SNP, falling back to 0.25."""
        # Try explicit MAF column first
        for col in ["maf", "alt_freq", "allele_frequency", "MAF"]:
            if col in snp.index and not pd.isna(snp.get(col)):
                try:
                    val = float(snp[col])
                    if 0 < val < 1:
                        return val
                except (ValueError, TypeError):
                    pass
        # Default: MAF=0.25 → E[G²] = 2*0.25*0.75 + 4*0.0625 = 0.625
        return 0.25


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Uncertainty Propagation Engine — 3-layer PRS uncertainty"
    )
    parser.add_argument("--prs-data", required=True,
                       help="PRS data CSV (raw or adjusted)")
    parser.add_argument("--ancestry", required=True,
                       help="Ancestry inference JSON")
    parser.add_argument("--snp-db", required=True,
                       help="Position-annotated SNP database CSV")
    parser.add_argument("--vcf", help="VCF for genotype uncertainty extraction")
    parser.add_argument("--consistency-check", help="Consistency check JSON")
    parser.add_argument("--output-dir", "-o", default="prs",
                       help="Output directory")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    prs_data = pd.read_csv(args.prs_data)

    engine = UncertaintyPropagationEngine(vcf_path=args.vcf)
    result = engine.propagate(
        prs_data=prs_data,
        ancestry_path=args.ancestry,
        snp_database_path=args.snp_db,
        output_dir=args.output_dir,
        consistency_check_path=args.consistency_check,
    )

    df = result["prs_uncertainty"]
    print("\n═══ Uncertainty-Quantified PRS ═══")
    print(f"Global uncertainty score: {result['report'].global_uncertainty_score:.4f}")
    print(f"\n{'Trait':<30} {'PRS':>8} {'±SE':>8} {'95% CI':>20} {'U-score':>8}")
    print("-" * 80)
    for _, row in df.iterrows():
        print(f"{row['trait']:<30} {row['prs']:>8.4f} {row['prs_se']:>8.4f} "
              f"[{row['ci_95_lower']:>8.4f}, {row['ci_95_upper']:<8.4f}] "
              f"{row['uncertainty_score']:>8.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
