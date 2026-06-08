#!/usr/bin/env python3
"""
PCA-Based Ancestry Classifier

Classifies target sample ancestry using 1000 Genomes PCA projection.
Uses genome-wide 2M LD-pruned variants (not 109 trait SNPs).

Methods:
  1. k-NN: majority vote of k nearest 1000G neighbors in PCA space
  2. Centroid: nearest population centroid (Mahalanobis distance)
  3. Ensemble: weighted vote of both methods

Input:
  pca/1000G_pcs.eigenvec     — 2,504 reference samples × 20 PCs
  pca/target_pcs.eigenvec    — 1 target sample × 20 PCs
  reference/1000G_full/population_panel.txt

Output:
  pca/ancestry_classification.json
  science/ANCESTRY_MODEL.json  (overwrites UNKNOWN with real data)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _venv

import json, logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SUPER_POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
POP_NAMES = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
             "SAS": "South Asian", "AMR": "Admixed American"}


def load_data(ref_pcs_path: str, target_pcs_path: str, pop_panel_path: str):
    """Load reference PCs, target PCs, and population labels."""
    ref = pd.read_csv(ref_pcs_path, sep=r"\s+")
    target = pd.read_csv(target_pcs_path, sep=r"\s+")
    panel = pd.read_csv(pop_panel_path, sep=r"\s+", dtype=str)

    pc_cols = [c for c in ref.columns if c.startswith("PC")]
    n_pcs = len(pc_cols)
    logger.info(f"  PCs available: {n_pcs}")

    # Map IID → super population
    pop_map = {}
    for _, row in panel.iterrows():
        pop_map[str(row.iloc[0])] = str(row.iloc[2])  # sample → super_pop

    ref["population"] = ref["IID"].astype(str).map(pop_map)
    ref = ref[ref["population"].isin(SUPER_POPS)]
    logger.info(f"  Reference: {len(ref)} samples, {ref['population'].nunique()} populations")

    # Extract PC matrices
    X_ref = ref[pc_cols].values.astype(np.float64)
    x_target = target[pc_cols].values[0].astype(np.float64)
    y_ref = ref["population"].values

    return X_ref, x_target, y_ref, pc_cols


def classify_knn(X_ref: np.ndarray, x_target: np.ndarray, y_ref: np.ndarray,
                 k: int = 15) -> Dict:
    """k-NN classification in PCA space."""
    from scipy.spatial import cKDTree
    tree = cKDTree(X_ref)
    distances, indices = tree.query(x_target.reshape(1, -1), k=k)

    # Weighted vote: 1/distance
    votes = {}
    for d, idx in zip(distances[0], indices[0]):
        pop = y_ref[idx]
        weight = 1.0 / max(d, 0.001)
        votes[pop] = votes.get(pop, 0) + weight

    total = sum(votes.values())
    probs = {pop: round(votes.get(pop, 0) / total, 4) for pop in SUPER_POPS}
    assigned = max(probs, key=probs.get)

    return {"method": "k-NN", "k": k, "assigned": assigned,
            "probabilities": probs, "mean_distance": float(np.mean(distances))}


def classify_centroid(X_ref: np.ndarray, x_target: np.ndarray, y_ref: np.ndarray) -> Dict:
    """Centroid distance classification."""
    centroids = {}
    for pop in SUPER_POPS:
        mask = y_ref == pop
        if mask.sum() > 0:
            centroids[pop] = X_ref[mask].mean(axis=0)

    distances = {}
    for pop, centroid in centroids.items():
        distances[pop] = float(np.linalg.norm(x_target - centroid))

    # Convert distances to probabilities (softmax of negative distances)
    max_dist = max(distances.values()) + 1.0
    probs = {}
    for pop in SUPER_POPS:
        d = distances.get(pop, max_dist)
        probs[pop] = round(np.exp(-d / max_dist), 4)

    total = sum(probs.values())
    probs = {pop: round(v / total, 4) for pop, v in probs.items()}
    assigned = max(probs, key=probs.get)

    return {"method": "Centroid", "assigned": assigned,
            "probabilities": probs, "distances": distances}


def classify_ensemble(knn_result: Dict, centroid_result: Dict) -> Dict:
    """Ensemble: average of k-NN and centroid probabilities."""
    probs = {}
    for pop in SUPER_POPS:
        probs[pop] = round(
            (knn_result["probabilities"].get(pop, 0) +
             centroid_result["probabilities"].get(pop, 0)) / 2, 4
        )
    assigned = max(probs, key=probs.get)

    # Confidence: top probability / second-best
    sorted_probs = sorted(probs.values(), reverse=True)
    top = sorted_probs[0]
    second = sorted_probs[1] if len(sorted_probs) > 1 else 0
    confidence_ratio = top / max(second, 0.01)

    if confidence_ratio > 5:
        confidence = "HIGH"
    elif confidence_ratio > 2:
        confidence = "MODERATE"
    else:
        confidence = "LOW"

    return {"method": "PCA_ENSEMBLE", "assigned": assigned,
            "probabilities": probs, "confidence": confidence,
            "confidence_ratio": round(confidence_ratio, 2)}


def classify(ref_pcs_path: str = "pca/1000G_pcs.eigenvec",
             target_pcs_path: str = "pca/target_pcs.eigenvec",
             pop_panel_path: str = "reference/1000G_full/population_panel.txt",
             output_dir: str = "pca") -> Dict:
    logger.info("═══ PCA-Based Ancestry Classification ═══")

    X_ref, x_target, y_ref, pc_cols = load_data(
        ref_pcs_path, target_pcs_path, pop_panel_path)

    # Run classifiers
    knn = classify_knn(X_ref, x_target, y_ref)
    centroid = classify_centroid(X_ref, x_target, y_ref)
    ensemble = classify_ensemble(knn, centroid)

    # Show results
    for pop in SUPER_POPS:
        prob = ensemble["probabilities"].get(pop, 0)
        bar = "█" * int(prob * 50)
        logger.info(f"    {pop}: {prob:.2%} {bar}")

    logger.info(f"  Assigned: {ensemble['assigned']} ({POP_NAMES.get(ensemble['assigned'], '')})")
    logger.info(f"  Confidence: {ensemble['confidence']} ({ensemble['confidence_ratio']:.1f}x)")

    # Save classification
    output_dir = Path(output_dir)
    result = {
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        "method": "PCA_ENSEMBLE_V2",
        "reference": "1000 Genomes Phase 3 (genome-wide, 2M LD-pruned variants, 20 PCs)",
        "n_pcs": len(pc_cols),
        "n_reference_samples": len(X_ref),
        "assigned_population": ensemble["assigned"],
        "assigned_population_name": POP_NAMES.get(ensemble["assigned"], ""),
        "confidence": ensemble["confidence"],
        "confidence_ratio": ensemble["confidence_ratio"],
        "posterior_probabilities": ensemble["probabilities"],
        "classifiers": {
            "knn": {"assigned": knn["assigned"], "probabilities": knn["probabilities"],
                    "mean_distance": knn["mean_distance"]},
            "centroid": {"assigned": centroid["assigned"], "probabilities": centroid["probabilities"]},
        },
        "is_valid_for_scoring": ensemble["confidence"] != "LOW",
    }

    with open(output_dir / "ancestry_classification.json", "w") as fh:
        json.dump(result, fh, indent=2)

    # Also update ANCESTRY_MODEL.json
    science_dir = Path("science")
    science_dir.mkdir(exist_ok=True)
    ssst_result = {
        "method": "PCA_ENSEMBLE_V2",
        "reference_panel": "1000 Genomes Phase 3 (all autosomes)",
        "n_pcs": len(pc_cols),
        "n_reference_samples": len(X_ref),
        "super_populations": SUPER_POPS,
        "assigned_population": ensemble["assigned"],
        "posterior_probabilities": ensemble["probabilities"],
        "confidence": ensemble["confidence"],
        "quality_metrics": {
            "confidence_ratio": ensemble["confidence_ratio"],
            "knn_mean_distance": knn["mean_distance"],
        },
        "is_valid_for_scoring": ensemble["confidence"] != "LOW",
        "model_hash": hex(hash(frozenset(ensemble["probabilities"].items())))[:16],
        "frozen_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
    }
    with open(science_dir / "ANCESTRY_MODEL.json", "w") as fh:
        json.dump(ssst_result, fh, indent=2)

    logger.info(f"  ✅ Classifier: {output_dir}/ancestry_classification.json")
    logger.info(f"  ✅ ANCESTRY_MODEL: science/ANCESTRY_MODEL.json")

    return result


def main():
    import argparse
    p = argparse.ArgumentParser(description="PCA-Based Ancestry Classifier")
    p.add_argument("--ref-pcs", default="pca/1000G_pcs.eigenvec")
    p.add_argument("--target-pcs", default="pca/target_pcs.eigenvec")
    p.add_argument("--pop-panel", default="reference/1000G_full/population_panel.txt")
    p.add_argument("--output-dir", "-o", default="pca")
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")
    return 0 if classify(args.ref_pcs, args.target_pcs, args.pop_panel, args.output_dir) else 1


if __name__ == "__main__":
    sys.exit(main())
