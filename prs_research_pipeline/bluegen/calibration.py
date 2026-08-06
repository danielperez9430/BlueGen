"""
Stage H — Population Calibration (IMPROVEMENT_PLAN.md 2.1).

Extracted from scripts/prs/population_calibrate_v2.py, which is now a thin
CLI wrapper around this module. Logic unchanged - see that file's git
history for the original bug fixes this module inherits (in particular the
n_samples<10 low_confidence fix documented in calibrate_sample() below).

For each trait x population:
  - Compute PRS for all 1000G samples in that population
  - Build empirical distribution: mean, sd, median, IQR, percentiles
  - Store as reference distribution

For target sample:
  - Population-specific z-score: z_pop = (PRS - mu_pop) / sigma_pop
  - Population-specific percentile: from empirical CDF
  - Global z-score: using pooled distribution

No hardcoded mu=0. No synthetic assumptions.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logger = logging.getLogger(__name__)

SUPER_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]


@dataclass
class PopulationDistribution:
    """Empirical PRS distribution for one trait × population."""
    trait: str
    population: str
    n_samples: int
    mean: float
    std: float
    median: float
    iqr: float
    percentile_5: float
    percentile_25: float
    percentile_75: float
    percentile_95: float
    skewness: float
    kurtosis: float
    shapiro_p: float  # Normality test p-value


@dataclass
class CalibratedPRS:
    """Population-calibrated PRS for a single trait."""
    trait: str
    prs_raw: float
    assigned_population: str
    population_mu: float
    population_sigma: float
    z_score_population: float
    percentile_population: float
    z_score_global: float
    percentile_global: float
    risk_category: str   # low/medium/high based on population percentiles
    low_confidence: bool = False   # True if the reference distribution is missing/too small (<10 samples)
    n_reference_samples: int = 0   # 1000G reference samples backing population_mu/population_sigma


class PopulationCalibrationV2:
    """
    Computes and applies empirical population-specific PRS reference distributions.

    Key correction vs. v1:
      - μ_pop computed from actual 1000G data, not hardcoded to 0
      - σ_pop computed from actual 1000G data, not hardcoded estimates
      - Empirical percentiles, not KDE on Gaussian simulations
      - Per-population sample sizes reported for transparency

    Usage:
        calibrator = PopulationCalibrationV2()
        calibrator.build_reference_distributions(
            ref_prs="prs/1000G_prs.csv",
            population_panel="reference/1000G_full/population_panel.txt",
            output_dir="reference/population_distributions/",
        )
        calibrated = calibrator.calibrate_sample(
            sample_prs="prs/prs_adjusted.csv",
            ancestry_json="ancestry/posterior_probabilities.json",
            output_dir="prs/",
        )
    """

    RISK_THRESHOLDS = {
        "low": 25,     # Bottom 25% of population
        "high": 75,    # Top 25% of population
    }

    def __init__(self):
        self._reference_distributions: Dict[str, Dict[str, PopulationDistribution]] = {}

    # ── Public API ───────────────────────────────────────────────────────

    def build_reference_distributions(
        self,
        ref_prs: str,
        population_panel: str,
        output_dir: str,
    ) -> Dict[str, Dict[str, PopulationDistribution]]:
        """
        Build empirical reference distributions from 1000 Genomes PRS data.

        Args:
            ref_prs: CSV with PRS computed for all 1000G samples.
                     Columns: individual_id, trait, prs_raw
            population_panel: 1000G population labels.
            output_dir: Output directory.

        Returns:
            Nested dict[trait][population] → PopulationDistribution
        """
        logger.info("═══ Building Empirical Population Distributions (Phase 6) ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load data
        prs_df = pd.read_csv(ref_prs)
        panel = pd.read_csv(population_panel, sep=r"\s+", dtype=str)

        # Build sample → population mapping
        sample_to_pop = {}
        for _, row in panel.iterrows():
            sample_to_pop[str(row.iloc[0])] = str(row.iloc[2]) if len(row.columns) >= 3 else str(row.iloc[1])

        # Map PRS samples to populations
        prs_df["population"] = prs_df["individual_id"].astype(str).map(sample_to_pop)
        prs_df = prs_df[prs_df["population"].isin(SUPER_POPULATIONS)]

        logger.info(f"  Reference PRS: {len(prs_df)} trait-samples across {len(SUPER_POPULATIONS)} populations")

        # Compute distributions per trait × population
        traits = sorted(prs_df["trait"].unique())
        distributions = {}

        for trait in traits:
            trait_data = prs_df[prs_df["trait"] == trait]
            distributions[trait] = {}

            for pop in SUPER_POPULATIONS:
                pop_data = trait_data[trait_data["population"] == pop]["prs_raw"]
                n = len(pop_data)

                if n < 5:
                    logger.warning(f"    {trait}/{pop}: insufficient samples (n={n})")
                    continue

                values = pop_data.values.astype(np.float64)

                # Compute distribution statistics
                mean = float(np.mean(values))
                std = float(np.std(values, ddof=1))
                median = float(np.median(values))
                q25, q75 = np.percentile(values, [25, 75])
                p5, p95 = np.percentile(values, [5, 95])
                skewness = float(scipy_stats.skew(values))
                kurtosis = float(scipy_stats.kurtosis(values))

                # Normality test
                if n >= 8:
                    _, shapiro_p = scipy_stats.shapiro(values[:min(n, 5000)])
                else:
                    shapiro_p = np.nan

                dist = PopulationDistribution(
                    trait=trait,
                    population=pop,
                    n_samples=n,
                    mean=round(mean, 6),
                    std=round(std, 6),
                    median=round(median, 6),
                    iqr=round(float(q75 - q25), 6),
                    percentile_5=round(float(p5), 6),
                    percentile_25=round(float(q25), 6),
                    percentile_75=round(float(q75), 6),
                    percentile_95=round(float(p95), 6),
                    skewness=round(skewness, 4),
                    kurtosis=round(kurtosis, 4),
                    shapiro_p=round(float(shapiro_p), 6),
                )

                distributions[trait][pop] = dist
                logger.info(f"    {trait}/{pop}: μ={mean:.3f}, σ={std:.3f}, n={n}")

        self._reference_distributions = distributions

        # Save
        self._save_distributions(distributions, output_dir)

        return distributions

    def calibrate_sample(
        self,
        sample_prs: str,
        ancestry_json: str,
        output_dir: str,
        sample_id: str = "SAMPLE_001",
        ref_dist_dir: Optional[str] = None,
        distributions: Optional[Dict[str, Dict[str, PopulationDistribution]]] = None,
    ) -> List[CalibratedPRS]:
        """
        Calibrate a sample's PRS against empirical population distributions.

        Args:
            sample_prs: Sample PRS CSV.
            ancestry_json: Ancestry classification JSON.
            output_dir: Output directory.
            sample_id: Sample identifier.
            ref_dist_dir: Directory with pre-computed reference distributions,
                loaded via _load_distributions() if `distributions` isn't
                passed directly. Kept for the CLI wrapper's build/load flow.
            distributions: Reference distributions to calibrate against,
                passed explicitly instead of relying on a prior
                build_reference_distributions()/_load_distributions() call
                having populated instance state - this is what makes
                calibration a pure function of its inputs, testable without
                reaching into a private attribute (IMPROVEMENT_PLAN.md 2.1).

        Returns:
            List of CalibratedPRS per trait.
        """
        logger.info("═══ Calibrating Sample PRS (Phase 6) ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Reference distributions: explicit param wins, else fall back to the
        # build/load-populated instance state (CLI wrapper's flow).
        if distributions is not None:
            self._reference_distributions = distributions
        elif ref_dist_dir:
            self._load_distributions(ref_dist_dir)

        if not self._reference_distributions:
            logger.error("  No reference distributions available. Run build_reference_distributions() first.")
            return []

        # Load ancestry
        with open(ancestry_json) as fh:
            ancestry = json.load(fh)
        assigned_pop = ancestry.get("assigned_population", "EUR")
        ancestry_probs = ancestry.get("posterior_probabilities", {assigned_pop: 1.0})

        # Load sample PRS
        prs_df = pd.read_csv(sample_prs)

        calibrated = []
        for _, row in prs_df.iterrows():
            trait = str(row.get("trait", ""))
            prs_raw = float(row.get("prs_adjusted", row.get("prs_raw", 0)))

            # Get population distribution
            pop_dist = self._get_distribution(trait, assigned_pop)

            if pop_dist is None:
                # No reference distribution at all: mu/sigma are meaningless,
                # so z and percentile must both stay neutral (not just one of
                # them - see IMPROVEMENT_PLAN.md TIER 3, calibration bug fix).
                pop_mu = 0.0
                pop_sigma = 1.0
                n_ref = 0
                low_confidence = True
                logger.warning(f"    No reference distribution for {trait}/{assigned_pop}")
                z_pop = 0.0
                percentile_pop = 50.0
            else:
                pop_mu = pop_dist.mean
                pop_sigma = pop_dist.std if pop_dist.std > 0 else 1.0
                n_ref = pop_dist.n_samples
                low_confidence = pop_dist.n_samples < 10

                # Population-specific z-score
                z_pop = (prs_raw - pop_mu) / pop_sigma if pop_sigma > 0 else 0.0

                # Population-specific percentile - always derived from z_pop
                # (via the normal approximation) so the two never disagree.
                # Previously this was hardcoded to 50.0 when n_samples < 10
                # while z_pop was still computed normally, which silently
                # forced risk_category to "medium" for those traits and fed
                # a degenerate calibration slope into 27_real_world_calibration.py.
                # low_confidence flags the same condition instead, so callers
                # can add an uncertainty caveat without discarding the signal.
                percentile_pop = scipy_stats.norm.cdf(z_pop) * 100

            # Global z-score (using all-population pooled stats)
            global_stats = self._compute_global_stats(trait, ancestry_probs)
            z_global = (prs_raw - global_stats["mu"]) / max(global_stats["sigma"], 0.001)
            percentile_global = scipy_stats.norm.cdf(z_global) * 100

            # Risk category
            if percentile_pop >= self.RISK_THRESHOLDS["high"]:
                risk = "high"
            elif percentile_pop <= self.RISK_THRESHOLDS["low"]:
                risk = "low"
            else:
                risk = "medium"

            calibrated.append(CalibratedPRS(
                trait=trait,
                prs_raw=round(prs_raw, 4),
                assigned_population=assigned_pop,
                population_mu=round(pop_mu, 4),
                population_sigma=round(pop_sigma, 4),
                z_score_population=round(z_pop, 4),
                percentile_population=round(percentile_pop, 1),
                z_score_global=round(z_global, 4),
                percentile_global=round(percentile_global, 1),
                risk_category=risk,
                low_confidence=low_confidence,
                n_reference_samples=n_ref,
            ))

        # Save
        self._save_calibrated(calibrated, sample_id, output_dir)

        return calibrated

    # ── Private: Distribution Management ──────────────────────────────────

    def _get_distribution(
        self, trait: str, population: str
    ) -> Optional[PopulationDistribution]:
        """Get reference distribution for a trait × population."""
        return self._reference_distributions.get(trait, {}).get(population)

    def _compute_global_stats(
        self, trait: str, ancestry_probs: Dict[str, float]
    ) -> Dict[str, float]:
        """Compute ancestry-weighted global mean and variance."""
        weighted_mu = 0.0
        weighted_var = 0.0
        total_weight = 0.0

        for pop, prob in ancestry_probs.items():
            dist = self._get_distribution(trait, pop)
            if dist is not None:
                weighted_mu += prob * dist.mean
                weighted_var += prob * (dist.std ** 2)
                total_weight += prob

        if total_weight > 0:
            weighted_mu /= total_weight
            weighted_var /= total_weight

        return {"mu": weighted_mu, "sigma": np.sqrt(max(weighted_var, 1e-6))}

    def _save_distributions(
        self,
        distributions: Dict[str, Dict[str, PopulationDistribution]],
        output_dir: Path,
    ) -> None:
        """Save reference distributions to disk."""
        # Save as JSON
        json_data = {}
        for trait, pops in distributions.items():
            json_data[trait] = {}
            for pop, dist in pops.items():
                json_data[trait][pop] = asdict(dist)

        json_path = output_dir / "reference_distributions.json"
        with open(json_path, "w") as fh:
            json.dump(json_data, fh, indent=2)

        # Save as CSV
        rows = []
        for trait, pops in distributions.items():
            for pop, dist in pops.items():
                rows.append(asdict(dist))
        csv_path = output_dir / "reference_distributions.csv"
        pd.DataFrame(rows).to_csv(csv_path, index=False)

        logger.info(f"  ✅ Distributions: {json_path} ({len(rows)} trait×population entries)")

    def _load_distributions(self, dist_dir: str) -> None:
        """Load pre-computed reference distributions."""
        json_path = Path(dist_dir) / "reference_distributions.json"
        if not json_path.exists():
            logger.warning(f"    Distributions not found: {json_path}")
            return

        with open(json_path) as fh:
            raw = json.load(fh)

        # Support both nested format (distributions inside a wrapper) and flat format
        dists = raw.get("distributions", raw)

        self._reference_distributions = {}
        for trait, pops in dists.items():
            if not isinstance(pops, dict):
                continue  # Skip metadata keys like "generated_date", "n_populations"
            self._reference_distributions[trait] = {}
            for pop, data in pops.items():
                if not isinstance(data, dict):
                    continue
                # Normalize alternate key names (q25→percentile_25, etc.)
                normalized = dict(data)
                key_map = {"q25": "percentile_25", "q75": "percentile_75",
                           "p5": "percentile_5", "p95": "percentile_95"}
                for old_key, new_key in key_map.items():
                    if old_key in normalized and new_key not in normalized:
                        normalized[new_key] = normalized.pop(old_key)
                # Compute iqr from percentile_25/75 if missing
                if "iqr" not in normalized and "percentile_25" in normalized and "percentile_75" in normalized:
                    normalized["iqr"] = normalized["percentile_75"] - normalized["percentile_25"]
                # Filter to only fields the dataclass accepts
                valid_fields = {f.name for f in PopulationDistribution.__dataclass_fields__.values()}
                normalized = {k: v for k, v in normalized.items() if k in valid_fields}
                self._reference_distributions[trait][pop] = PopulationDistribution(**normalized)

        logger.info(f"  Loaded distributions for {len(self._reference_distributions)} traits")

    def _save_calibrated(
        self, calibrated: List[CalibratedPRS], sample_id: str, output_dir: Path
    ) -> None:
        """Save calibrated PRS to CSV."""
        rows = []
        for c in calibrated:
            rows.append({
                "individual_id": sample_id,
                "trait": c.trait,
                "prs_raw": c.prs_raw,
                "assigned_population": c.assigned_population,
                "population_mu": c.population_mu,
                "population_sigma": c.population_sigma,
                "z_score_population": c.z_score_population,
                "percentile_population": c.percentile_population,
                "z_score_global": c.z_score_global,
                "percentile_global": c.percentile_global,
                "risk_category": c.risk_category,
                "low_confidence": c.low_confidence,
                "n_reference_samples": c.n_reference_samples,
            })

        df = pd.DataFrame(rows)
        csv_path = output_dir / "population_calibrated_v2.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"  ✅ Calibrated PRS: {csv_path}")

        # Log summary
        for c in calibrated:
            pop_flag = f" [{c.assigned_population}]" if c.assigned_population != "EUR" else ""
            logger.info(f"    {c.trait}: z={c.z_score_population:+.2f}{pop_flag}, "
                      f"pctl={c.percentile_population:.0f}%, risk={c.risk_category}")
