#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 6 — MODULE 3: REAL ANCESTRY INFERENCE                               ║
║   scripts/ancestry_inference_v2.py                                          ║
║                                                                            ║
║   Replaces the 33-SNP allele-frequency method with proper PCA-based         ║
║   ancestry classification using genome-wide markers.                        ║
║                                                                            ║
║   Methods (all using genome-wide PCA coordinates):                          ║
║     1. Centroid distance — Euclidean distance to population centroids       ║
║     2. Multinomial logistic regression — trained on 1000G labels            ║
║     3. Nearest-neighbor consensus — majority vote among k nearest refs     ║
║                                                                            ║
║   Super-populations: EUR, AFR, EAS, SAS, AMR                                ║
║                                                                            ║
║   Quality metrics:                                                          ║
║     • Posterior probability (softmax over distances)                        ║
║     • Classification entropy (uncertainty)                                  ║
║     • Nearest population distance                                           ║
║     • Second-best distance ratio (discrimination)                           ║
║     • Confidence flag (HIGH/MODERATE/LOW/REJECT)                            ║
║                                                                            ║
║   Rejects ancestry calls with insufficient evidence.                        ║
║                                                                            ║
║   Key correction: No longer uses trait-associated SNPs for ancestry.        ║
║                                                                            ║
║   Output:                                                                   ║
║     ancestry/posterior_probabilities.json                                   ║
║     ancestry/classification_report.json                                     ║
║     ancestry/quality_metrics.json                                           ║
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
from scipy.spatial.distance import cdist, mahalanobis

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

SUPER_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
POP_NAMES = {
    "EUR": {"en": "European", "es": "Europea"},
    "AFR": {"en": "African", "es": "Africana"},
    "EAS": {"en": "East Asian", "es": "Asiática Oriental"},
    "SAS": {"en": "South Asian", "es": "Sudasiática"},
    "AMR": {"en": "Admixed American", "es": "Americana Mixta"},
}

CONFIDENCE_THRESHOLDS = {
    "HIGH": 0.90,      # ≥90% posterior probability
    "MODERATE": 0.70,  # 70-90%
    "LOW": 0.50,       # 50-70%
    # <50% → REJECT
}


# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class AncestryClassification:
    """Complete ancestry classification for a sample."""
    sample_id: str
    assigned_population: str
    confidence: str                  # HIGH, MODERATE, LOW, REJECT
    posterior_probabilities: Dict[str, float]
    method_agreement: Dict[str, str]  # method → assigned population
    quality_metrics: Dict[str, float]
    centroid_distances: Dict[str, float]
    nearest_population: str
    second_nearest: str
    distance_ratio: float            # d2/d1 — discrimination quality
    entropy: float                   # Shannon entropy of probabilities


# ── Ancestry Inference Engine V2 ──────────────────────────────────────────────

class AncestryInferenceV2:
    """
    PCA-based ancestry inference using genome-wide markers.

    Three independent methods are combined:
      1. Centroid distance (Euclidean in PC space)
      2. Logistic regression (multinomial, trained on 1000G)
      3. k-Nearest Neighbors (k=50, majority vote)

    The ensemble provides:
      - Method agreement checks (discordant methods → lower confidence)
      - Posterior probabilities via softmax over centroid distances
      - Quality metrics for rejection of ambiguous calls

    Usage:
        inferrer = AncestryInferenceV2()
        result = inferrer.classify(
            projected_pcs="pca/projected_sample.csv",
            ref_pcs="pca/1000G_pca.eigenvec",
            ref_centroids="pca/reference_centroids.csv",
            population_panel="reference/1000G_full/population_panel.txt",
            output_dir="ancestry/",
        )
    """

    def __init__(self):
        pass

    # ── Public API ───────────────────────────────────────────────────────

    def classify(
        self,
        projected_pcs: str,
        ref_pcs: str,
        ref_centroids: str,
        population_panel: str,
        output_dir: str,
        sample_id: str = "SAMPLE_001",
        n_pcs: int = 10,
        k_neighbors: int = 50,
    ) -> AncestryClassification:
        """
        Classify sample ancestry using PCA-based methods.

        Args:
            projected_pcs: CSV with sample PC coordinates.
            ref_pcs: 1000G reference PC coordinates (eigenvec).
            ref_centroids: Population centroids CSV.
            population_panel: 1000G population labels.
            output_dir: Output directory.
            sample_id: Sample identifier.
            n_pcs: Number of PCs to use.
            k_neighbors: Number of neighbors for k-NN.

        Returns:
            AncestryClassification with probabilities and quality metrics.
        """
        logger.info("═══ Ancestry Inference V2 (Phase 6) ═══")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load data
        sample_pcs = self._load_sample_pcs(projected_pcs, n_pcs)
        ref_coords, ref_labels = self._load_reference(ref_pcs, population_panel, n_pcs)
        centroids = self._load_centroids(ref_centroids, n_pcs)

        logger.info(f"  Reference samples: {len(ref_coords)}")
        logger.info(f"  Populations: {sorted(set(ref_labels))}")

        # Method 1: Centroid distance
        centroid_probs, centroid_distances = self._centroid_classify(
            sample_pcs, centroids
        )
        centroid_pop = max(centroid_probs, key=centroid_probs.get)

        # Method 2: Logistic regression
        lr_probs, lr_pop = self._logistic_classify(
            sample_pcs, ref_coords, ref_labels
        )

        # Method 3: k-Nearest Neighbors
        knn_pop = self._knn_classify(
            sample_pcs, ref_coords, ref_labels, k=k_neighbors
        )

        # Ensemble: weighted average of probabilities
        ensemble_probs = {}
        for pop in SUPER_POPULATIONS:
            ensemble_probs[pop] = (
                0.5 * centroid_probs.get(pop, 0) +
                0.3 * lr_probs.get(pop, 0) +
                0.2 * (1.0 if knn_pop == pop else 0.0)
            )

        # Assign population
        assigned_pop = max(ensemble_probs, key=ensemble_probs.get)
        max_prob = ensemble_probs[assigned_pop]

        # Determine confidence
        confidence = "REJECT"
        for level, threshold in sorted(CONFIDENCE_THRESHOLDS.items(),
                                       key=lambda x: x[1], reverse=True):
            if max_prob >= threshold:
                confidence = level
                break

        # Method agreement
        method_agreement = {
            "centroid": centroid_pop,
            "logistic_regression": lr_pop,
            "knn": knn_pop,
        }
        methods_agree = len(set(method_agreement.values())) == 1

        # Quality metrics
        sorted_pops = sorted(ensemble_probs.items(), key=lambda x: x[1], reverse=True)
        second_nearest = sorted_pops[1][0] if len(sorted_pops) > 1 else "UNKNOWN"
        d1 = centroid_distances.get(assigned_pop, 1.0)
        d2 = centroid_distances.get(second_nearest, d1 * 2)
        distance_ratio = d2 / max(d1, 1e-10)

        # Shannon entropy of probabilities
        probs_array = np.array([ensemble_probs.get(p, 0.01) for p in SUPER_POPULATIONS])
        probs_array = probs_array / probs_array.sum()
        entropy = -np.sum(probs_array * np.log(np.maximum(probs_array, 1e-10)))

        quality_metrics = {
            "max_posterior_probability": round(float(max_prob), 4),
            "entropy": round(float(entropy), 4),
            "distance_ratio": round(float(distance_ratio), 4),
            "methods_agree": methods_agree,
            "n_methods_agreeing": len(set(method_agreement.values())),
            "n_reference_samples": len(ref_coords),
            "n_pcs_used": n_pcs,
        }

        # Build result
        result = AncestryClassification(
            sample_id=sample_id,
            assigned_population=assigned_pop,
            confidence=confidence,
            posterior_probabilities={k: round(v, 4) for k, v in ensemble_probs.items()},
            method_agreement=method_agreement,
            quality_metrics=quality_metrics,
            centroid_distances={k: round(v, 4) for k, v in centroid_distances.items()},
            nearest_population=assigned_pop,
            second_nearest=second_nearest,
            distance_ratio=round(float(distance_ratio), 4),
            entropy=round(float(entropy), 4),
        )

        # Log result
        agreement_str = "✅ ALL AGREE" if methods_agree else f"⚠️ DISCORDANT: {method_agreement}"
        logger.info(f"  Assigned: {assigned_pop} ({confidence}, p={max_prob:.3f})")
        logger.info(f"  Methods: {agreement_str}")
        logger.info(f"  Distance ratio: {distance_ratio:.2f} (d2/d1)")

        if confidence == "REJECT":
            logger.warning(f"  ⚠️  ANCESTRY CALL REJECTED — insufficient confidence ({max_prob:.1%})")

        # Save outputs
        self._save_outputs(result, output_dir)

        return result

    # ── Method 1: Centroid Distance ──────────────────────────────────────

    def _centroid_classify(
        self, sample_pcs: np.ndarray, centroids: Dict[str, np.ndarray]
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        """Classify by Euclidean distance to population centroids."""
        distances = {}
        for pop, centroid in centroids.items():
            diff = sample_pcs - centroid
            distances[pop] = float(np.sqrt(np.sum(diff ** 2)))

        if not distances:
            return {}, {}

        # Softmax over negative distances (closer = higher probability)
        dist_array = np.array([distances[p] for p in SUPER_POPULATIONS if p in distances])
        pops_present = [p for p in SUPER_POPULATIONS if p in distances]

        # Temperature scaling: smaller temperature = sharper distribution
        temperature = np.median(dist_array) * 0.5 if np.median(dist_array) > 0 else 1.0
        neg_dists = -dist_array / temperature
        exp_vals = np.exp(neg_dists - np.max(neg_dists))  # Stable softmax
        probs_array = exp_vals / exp_vals.sum()

        probs = {pop: float(probs_array[i]) for i, pop in enumerate(pops_present)}
        return probs, distances

    # ── Method 2: Logistic Regression ────────────────────────────────────

    def _logistic_classify(
        self,
        sample_pcs: np.ndarray,
        ref_coords: np.ndarray,
        ref_labels: np.ndarray,
    ) -> Tuple[Dict[str, float], str]:
        """Classify using multinomial logistic regression."""
        try:
            from sklearn.linear_model import LogisticRegression
        except ImportError:
            logger.warning("    scikit-learn not available — skipping logistic regression")
            return {}, "UNKNOWN"

        # Filter to populations with sufficient samples
        unique_labels, counts = np.unique(ref_labels, return_counts=True)
        valid_pops = unique_labels[counts >= 10]  # Need at least 10 samples per class

        valid_mask = np.isin(ref_labels, valid_pops)
        X_train = ref_coords[valid_mask]
        y_train = ref_labels[valid_mask]

        if len(valid_pops) < 2:
            return {}, str(valid_pops[0]) if len(valid_pops) == 1 else "UNKNOWN"

        # Train multinomial logistic regression
        clf = LogisticRegression(
            multi_class="multinomial",
            solver="lbfgs",
            max_iter=1000,
            C=1.0,  # L2 regularization
            class_weight="balanced",
            random_state=42,
        )
        clf.fit(X_train, y_train)

        # Predict probabilities
        sample_2d = sample_pcs.reshape(1, -1)
        proba = clf.predict_proba(sample_2d)[0]
        predicted = clf.predict(sample_2d)[0]

        probs = {pop: 0.0 for pop in SUPER_POPULATIONS}
        for i, pop in enumerate(clf.classes_):
            probs[str(pop)] = float(proba[i])

        return probs, str(predicted)

    # ── Method 3: k-Nearest Neighbors ────────────────────────────────────

    def _knn_classify(
        self,
        sample_pcs: np.ndarray,
        ref_coords: np.ndarray,
        ref_labels: np.ndarray,
        k: int = 50,
    ) -> str:
        """Classify by majority vote among k nearest reference neighbors."""
        # Compute distances to all reference samples
        distances = cdist(sample_pcs.reshape(1, -1), ref_coords, metric="euclidean")[0]

        # Get k nearest
        k = min(k, len(distances))
        nearest_idx = np.argpartition(distances, k)[:k]
        nearest_labels = ref_labels[nearest_idx]

        # Majority vote
        unique, counts = np.unique(nearest_labels, return_counts=True)
        winner = unique[np.argmax(counts)]

        return str(winner)

    # ── Data Loading ─────────────────────────────────────────────────────

    def _load_sample_pcs(self, path: str, n_pcs: int) -> np.ndarray:
        """Load sample PC coordinates from CSV."""
        df = pd.read_csv(path)
        pc_cols = [f"PC{i+1}" for i in range(n_pcs) if f"PC{i+1}" in df.columns]
        if not pc_cols:
            raise ValueError(f"No PC columns found in {path}")
        return df[pc_cols].values[0].astype(np.float64)

    def _load_reference(
        self, ref_pcs_path: str, panel_path: str, n_pcs: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Load reference PC coordinates and population labels."""
        # Load eigenvec
        eigenvec = pd.read_csv(ref_pcs_path, sep=r"\s+", header=None)
        pc_cols = list(range(2, min(2 + n_pcs, eigenvec.shape[1])))
        ref_coords = eigenvec.iloc[:, pc_cols].values.astype(np.float64)
        sample_ids = eigenvec.iloc[:, 1].astype(str).values

        # Load population panel
        panel = pd.read_csv(panel_path, sep=r"\s+", dtype=str)
        pop_map = {}
        for _, row in panel.iterrows():
            pop_map[str(row.iloc[0])] = str(row.iloc[2]) if len(row.columns) >= 3 else str(row.iloc[1])

        # Map labels
        ref_labels = np.array([pop_map.get(sid, "UNKNOWN") for sid in sample_ids])

        # Filter to known populations
        known_mask = np.isin(ref_labels, SUPER_POPULATIONS)
        logger.info(f"    Reference samples with known populations: {known_mask.sum()}/{len(ref_labels)}")

        return ref_coords[known_mask], ref_labels[known_mask]

    def _load_centroids(self, path: str, n_pcs: int) -> Dict[str, np.ndarray]:
        """Load population centroids from CSV."""
        df = pd.read_csv(path)
        centroids = {}
        pc_cols = [f"PC{i+1}" for i in range(n_pcs) if f"PC{i+1}" in df.columns]
        for _, row in df.iterrows():
            pop = str(row["population"])
            centroids[pop] = row[pc_cols].values.astype(np.float64)
        return centroids

    # ── Output ───────────────────────────────────────────────────────────

    def _save_outputs(
        self, result: AncestryClassification, output_dir: Path
    ) -> None:
        """Save all ancestry outputs."""

        # Posterior probabilities
        probs_path = output_dir / "posterior_probabilities.json"
        with open(probs_path, "w") as fh:
            json.dump({
                "sample_id": result.sample_id,
                "assigned_population": result.assigned_population,
                "confidence": result.confidence,
                "posterior_probabilities": result.posterior_probabilities,
                "method_agreement": result.method_agreement,
            }, fh, indent=2)
        logger.info(f"  ✅ Probabilities: {probs_path}")

        # Classification report
        report_path = output_dir / "classification_report.json"
        with open(report_path, "w") as fh:
            json.dump({
                "sample_id": result.sample_id,
                "classification": {
                    "assigned_population": result.assigned_population,
                    "confidence": result.confidence,
                    "population_name": POP_NAMES.get(result.assigned_population, {}),
                },
                "ensemble_probabilities": result.posterior_probabilities,
                "method_results": result.method_agreement,
                "methods_agree": result.quality_metrics.get("methods_agree", False),
                "recommendation": self._get_recommendation(result),
            }, fh, indent=2)
        logger.info(f"  ✅ Report: {report_path}")

        # Quality metrics
        quality_path = output_dir / "quality_metrics.json"
        with open(quality_path, "w") as fh:
            json.dump({
                "sample_id": result.sample_id,
                "quality_metrics": result.quality_metrics,
                "centroid_distances": result.centroid_distances,
                "distance_ratio": result.distance_ratio,
                "entropy": result.entropy,
                "validation": {
                    "confidence_thresholds": CONFIDENCE_THRESHOLDS,
                    "n_reference_populations": len(SUPER_POPULATIONS),
                    "reference_panel": "1000 Genomes Phase 3",
                    "method": "PCA_ensemble_v2",
                },
            }, fh, indent=2)
        logger.info(f"  ✅ Quality: {quality_path}")

    def _get_recommendation(self, result: AncestryClassification) -> Dict[str, str]:
        """Generate population-specific calibration recommendation."""
        if result.confidence == "REJECT":
            return {
                "en": "Ancestry classification rejected due to insufficient confidence. "
                      "Use global (uncalibrated) PRS percentiles only.",
                "es": "Clasificación de ascendencia rechazada por confianza insuficiente. "
                      "Usar solo percentiles PRS globales (no calibrados).",
            }

        pop = result.assigned_population
        pop_name = POP_NAMES.get(pop, {}).get("en", pop)

        if pop == "EUR":
            return {
                "en": f"Classified as {pop_name} ({result.confidence} confidence). "
                      "GWAS effect sizes are well-calibrated for European populations. "
                      "Population-specific reference distributions can be applied.",
                "es": f"Clasificado como {POP_NAMES.get(pop, {}).get('es', pop)} "
                      f"(confianza {result.confidence}). Los tamaños de efecto GWAS están "
                      "bien calibrados para poblaciones europeas.",
            }
        elif pop == "AFR":
            return {
                "en": f"Classified as {pop_name} ({result.confidence} confidence). "
                      "CAUTION: Most GWAS are European-derived. PRS transferability "
                      "to African populations is reduced. Use with explicit caveats.",
                "es": f"Clasificado como {POP_NAMES.get(pop, {}).get('es', pop)} "
                      f"(confianza {result.confidence}). PRECAUCIÓN: La mayoría de GWAS "
                      "son de origen europeo. Transferibilidad reducida a poblaciones africanas.",
            }
        else:
            return {
                "en": f"Classified as {pop_name} ({result.confidence} confidence). "
                      "GWAS transferability may be limited. Use population-calibrated "
                      "percentiles with documented caveats.",
                "es": f"Clasificado como {POP_NAMES.get(pop, {}).get('es', pop)} "
                      f"(confianza {result.confidence}). La transferibilidad GWAS puede "
                      "ser limitada. Usar percentiles calibrados con advertencias documentadas.",
            }


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 6 Module 3: Real Ancestry Inference (PCA-based, genome-wide)"
    )
    parser.add_argument("--projected-pcs", required=True,
                       help="Sample PC projection CSV from pca_true_projection.py")
    parser.add_argument("--ref-pcs", required=True,
                       help="1000G reference eigenvec file")
    parser.add_argument("--ref-centroids", required=True,
                       help="Population centroids CSV")
    parser.add_argument("--population-panel", required=True,
                       help="1000G population panel file")
    parser.add_argument("--output-dir", "-o", default="ancestry")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--n-pcs", type=int, default=10)
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    inferrer = AncestryInferenceV2()
    result = inferrer.classify(
        projected_pcs=args.projected_pcs,
        ref_pcs=args.ref_pcs,
        ref_centroids=args.ref_centroids,
        population_panel=args.population_panel,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        n_pcs=args.n_pcs,
    )

    print(f"\n═══ Ancestry Classification ═══")
    print(f"  Population: {result.assigned_population} ({result.confidence})")
    print(f"  Probabilities: {json.dumps(result.posterior_probabilities, indent=2)}")
    print(f"  Distance ratio: {result.distance_ratio:.2f}")
    print(f"  Entropy: {result.entropy:.4f}")
    print(f"  Methods: {result.method_agreement}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
