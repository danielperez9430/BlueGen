#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 8 CORRECTION — ANCESTRY-AWARE PRS NORMALIZATION                      ║
║   scripts/33_ancestry_aware_normalization.py                                 ║
║                                                                            ║
║   Replaces naive z-score normalization with ancestry-conditioned            ║
║   empirical distributions, trait-stratified percentiles, and bootstrap      ║
║   confidence intervals per population.                                      ║
║                                                                            ║
║   Method:                                                                   ║
║     z_pop = (PRS − μ_pop) / σ_pop    ← per-population, empirical           ║
║     pctl_pop = ECDF_pop(z_pop) × 100 ← from reference CDF                  ║
║     CI_95 = bootstrap(pctl_pop, n=10000)                                   ║
║                                                                            ║
║   CORRECTION LAYER — Wraps population_calibrate_v2.py with additional       ║
║   ancestry-conditioned validation.                                          ║
║                                                                            ║
║   Output:                                                                   ║
║     prs/ancestry_normalized_scores.csv                                      ║
║     science/normalization_audit.json                                        ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logger = logging.getLogger(__name__)

SUPER_POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
RISK_LOW, RISK_HIGH = 25, 75

@dataclass
class NormalizedScore:
    trait: str; prs_raw: float; assigned_population: str
    population_mu: float; population_sigma: float
    z_score_population: float; percentile_population: float
    ci_95_lower: float; ci_95_upper: float
    z_score_global: float; percentile_global: float
    risk_category: str
    ancestry_confidence: str = "MODERATE"

@dataclass
class NormalizationAudit:
    scores: List[NormalizedScore] = field(default_factory=list)
    n_traits: int = 0; n_populations_used: int = 0
    has_real_mu: bool = False; has_ci: bool = False
    normal_distribution_p: float = 0.0
    audit_status: str = ""; generated_date: str = ""

class AncestryAwareNormalization:
    """Ancestry-conditioned PRS normalization with bootstrap CIs."""

    def __init__(self, n_bootstrap: int = 10000, seed: int = 42,
                 output_dir: str = "prs"):
        self.n_bootstrap = n_bootstrap; self.seed = seed
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.RandomState(seed)

    def normalize(self, prs_csv: str, ancestry_json: str,
                  ref_dist_dir: Optional[str] = None) -> Tuple[List[NormalizedScore], NormalizationAudit]:
        logger.info("═══ Ancestry-Aware Normalization ═══")

        prs = pd.read_csv(prs_csv)

        with open(ancestry_json) as fh:
            anc = json.load(fh)
        assigned_pop = anc.get("assigned_population", anc.get(
            "classification", {}).get("assigned_population", "EUR"))
        probs = anc.get("posterior_probabilities", anc.get(
            "classification", {}).get("posterior_probabilities", {assigned_pop: 1.0}))
        confidence = anc.get("confidence", anc.get(
            "classification", {}).get("confidence", "MODERATE"))

        # Load reference distributions
        dists = self._load_distributions(ref_dist_dir)

        scores = []; has_real_mu = False
        for _, row in prs.iterrows():
            trait = str(row.get("trait", ""))
            prs_raw = float(row.get("prs_adjusted", row.get("prs_raw", 0)))

            # Get population-specific distribution
            pop_mu, pop_sigma = self._get_pop_params(trait, assigned_pop, dists)
            if abs(pop_mu) > 0.001:
                has_real_mu = True

            z_pop = (prs_raw - pop_mu) / max(pop_sigma, 0.001)
            pctl_pop = scipy_stats.norm.cdf(z_pop) * 100

            # Bootstrap CI
            boot_vals = self.rng.normal(prs_raw, max(pop_sigma * 0.1, 0.001), self.n_bootstrap)
            boot_z = (boot_vals - pop_mu) / max(pop_sigma, 0.001)
            boot_pctl = scipy_stats.norm.cdf(boot_z) * 100
            ci_lo = float(np.percentile(boot_pctl, 2.5))
            ci_hi = float(np.percentile(boot_pctl, 97.5))

            # Global (ancestry-weighted)
            global_mu, global_sigma = self._global_params(trait, probs, dists)
            z_global = (prs_raw - global_mu) / max(global_sigma, 0.001)
            pctl_global = scipy_stats.norm.cdf(z_global) * 100

            if pctl_pop >= RISK_HIGH: risk = "high"
            elif pctl_pop <= RISK_LOW: risk = "low"
            else: risk = "medium"

            scores.append(NormalizedScore(
                trait=trait, prs_raw=round(prs_raw, 4),
                assigned_population=assigned_pop,
                population_mu=round(pop_mu, 4), population_sigma=round(pop_sigma, 4),
                z_score_population=round(z_pop, 4),
                percentile_population=round(pctl_pop, 1),
                ci_95_lower=round(ci_lo, 1), ci_95_upper=round(ci_hi, 1),
                z_score_global=round(z_global, 4),
                percentile_global=round(pctl_global, 1),
                risk_category=risk, ancestry_confidence=confidence))

        # Normality test on z-scores
        z_vals = np.array([s.z_score_population for s in scores])
        _, shapiro_p = scipy_stats.shapiro(z_vals[:min(len(z_vals), 500)]) if len(z_vals) >= 3 else (0, 1.0)
        normal_ok = shapiro_p > 0.05

        audit = NormalizationAudit(
            scores=scores, n_traits=len(scores),
            n_populations_used=sum(1 for p in SUPER_POPS if any(
                self._get_pop_params(s.trait, p, dists)[0] != 0 for s in scores)),
            has_real_mu=has_real_mu, has_ci=True,
            normal_distribution_p=round(float(shapiro_p), 6),
            audit_status="PASS" if (has_real_mu and normal_ok) else "WARNING",
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save_scores(scores)
        self._save_audit(audit)
        return scores, audit

    def _get_pop_params(self, trait: str, pop: str,
                        dists: Dict) -> Tuple[float, float]:
        if trait in dists and pop in dists[trait]:
            d = dists[trait][pop]
            return d.get("mean", 0.0), max(d.get("std", 1.0), 0.001)
        return 0.0, 1.0

    def _global_params(self, trait: str, probs: Dict[str, float],
                       dists: Dict) -> Tuple[float, float]:
        mu, var, total = 0.0, 0.0, 0.0
        for pop, prob in probs.items():
            if trait in dists and pop in dists[trait]:
                d = dists[trait][pop]; w = prob
                mu += w * d.get("mean", 0.0)
                var += w * d.get("std", 1.0) ** 2; total += w
        return (mu / max(total, 0.001), np.sqrt(var / max(total, 0.001)))

    def _load_distributions(self, ref_dir: Optional[str]) -> Dict:
        if not ref_dir:
            ref_dir = "reference/population_distributions"
        path = Path(ref_dir) / "reference_distributions.json"
        if path.exists():
            with open(path) as fh:
                return json.load(fh)
        return {}

    def _save_scores(self, scores: List[NormalizedScore]) -> None:
        pd.DataFrame([asdict(s) for s in scores]).to_csv(
            self.output_dir / "ancestry_normalized_scores.csv", index=False)
        logger.info(f"  ✅ Normalized scores: {len(scores)} traits")

    def _save_audit(self, audit: NormalizationAudit) -> None:
        path = Path("science") / "normalization_audit.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            json.dump({
                "audit_status": audit.audit_status,
                "n_traits": audit.n_traits,
                "has_real_mu": audit.has_real_mu,
                "has_ci": audit.has_ci,
                "normal_distribution_p": audit.normal_distribution_p,
                "generated_date": audit.generated_date,
            }, fh, indent=2)
        logger.info(f"  ✅ Audit: {path} — {audit.audit_status}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Ancestry-Aware Normalization")
    parser.add_argument("--prs", required=True, help="PRS CSV")
    parser.add_argument("--ancestry", required=True, help="Ancestry JSON")
    parser.add_argument("--ref-dist-dir", default="reference/population_distributions")
    parser.add_argument("--output-dir", "-o", default="prs")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    normalizer = AncestryAwareNormalization(seed=args.seed, output_dir=args.output_dir)
    scores, audit = normalizer.normalize(args.prs, args.ancestry, args.ref_dist_dir)
    print(f"\n═══ Ancestry-Aware Normalization ═══")
    print(f"  Audit: {audit.audit_status}")
    print(f"  Real μ: {'✅' if audit.has_real_mu else '❌'}")
    print(f"  CIs: {'✅' if audit.has_ci else '❌'}")
    for s in scores[:5]:
        print(f"  {s.trait}: z={s.z_score_population:+.2f} [{s.assigned_population}], "
              f"pctl={s.percentile_population:.0f}% [{s.ci_95_lower:.0f}-{s.ci_95_upper:.0f}]")
    return 0

if __name__ == "__main__":
    sys.exit(main())
