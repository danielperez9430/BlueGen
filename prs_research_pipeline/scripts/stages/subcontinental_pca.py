#!/usr/bin/env python3
"""
Sub-Continental Ancestry Classifier (EUR fine-structure)

Classifies a EUR-assigned sample into European sub-populations (IBS/GBR/CEU/TSI/FIN)
using existing 1000 Genomes PCA coordinates and sub-population labels.

Unlike the continental classifier (5 super-populations), this one works only on
the EUR subset of 1000G (503 samples) and classifies into the 5 EUR sub-populations
that 1000G provides: IBS (Iberian), GBR (British), CEU (NW European), TSI (Tuscan), FIN (Finnish).

This uses the ALREADY-COMPUTED PCA coordinates — no new PLINK run needed.
The 20 PCs from the continental PCA capture genome-wide ancestry; subsetting to
EUR-only reference samples and re-classifying is sufficient for sub-population resolution.

Input:
  pca/1000G_pcs.eigenvec     — 2,504 reference samples × 20 PCs
  pca/target_pcs.eigenvec    — 1 target sample × 20 PCs
  reference/1000G_full/population_panel.txt

Output:
  pca/subcontinental_assignment.json
"""

import sys
import os
import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

EUR_SUB_POPS = ["IBS", "GBR", "CEU", "TSI", "FIN"]
SUB_POP_NAMES = {
    "IBS": "Iberian (Spain/Portugal)",
    "GBR": "British (England/Scotland)",
    "CEU": "Northwest European (Utah/German)",
    "TSI": "Tuscan Italian",
    "FIN": "Finnish",
}
SUB_POP_DESCRIPTIONS = {
    "IBS": "Southwestern European. High genetic affinity to Basques, Sardinians, and North Africans.",
    "GBR": "Northwestern European. Mix of Celtic, Anglo-Saxon, and Viking ancestry.",
    "CEU": "Central/Western European. Representative of continental Germanic populations.",
    "TSI": "Southern European. Bridges European and Mediterranean genetic clusters.",
    "FIN": "Northeastern European. Genetically distinct due to founder effects and Uralic admixture.",
}


def load_data(ref_pcs_path: str, target_pcs_path: str, pop_panel_path: str):
    """Load reference and target PCA data with sub-population labels."""
    ref = pd.read_csv(ref_pcs_path, sep=r"\s+")
    target = pd.read_csv(target_pcs_path, sep=r"\s+")
    panel = pd.read_csv(pop_panel_path, sep=r"\s+", dtype=str)

    pc_cols = [c for c in ref.columns if c.startswith("PC")]
    n_pcs = len(pc_cols)
    logger.info(f"  PCs available: {n_pcs}")

    # Map IID → sub_population (pop column) and super_population
    pop_map = {}
    super_map = {}
    for _, row in panel.iterrows():
        sid = str(row.iloc[0])
        pop_map[sid] = str(row.iloc[1])      # pop (e.g., IBS)
        super_map[sid] = str(row.iloc[2])     # super_pop (e.g., EUR)

    ref["sub_pop"] = ref["IID"].astype(str).map(pop_map)
    ref["super_pop"] = ref["IID"].astype(str).map(super_map)

    # Filter to EUR only
    ref_eur = ref[ref["super_pop"] == "EUR"]
    ref_eur = ref_eur[ref_eur["sub_pop"].isin(EUR_SUB_POPS)]
    logger.info(f"  EUR reference: {len(ref_eur)} samples across {ref_eur['sub_pop'].nunique()} sub-populations")
    for pop in EUR_SUB_POPS:
        count = len(ref_eur[ref_eur["sub_pop"] == pop])
        if count > 0:
            logger.info(f"    {pop}: {count} samples")

    X_ref = ref_eur[pc_cols].values.astype(np.float64)
    x_target = target[pc_cols].values[0].astype(np.float64)
    y_ref = ref_eur["sub_pop"].values

    return X_ref, x_target, y_ref, pc_cols


def classify_centroid(X_ref: np.ndarray, x_target: np.ndarray, y_ref: np.ndarray,
                      pops: list) -> Dict:
    """Classify by nearest population centroid (Euclidean distance in PCA space)."""
    centroids = {}
    for pop in pops:
        mask = y_ref == pop
        if mask.sum() > 0:
            centroids[pop] = X_ref[mask].mean(axis=0)

    distances = {}
    for pop, centroid in centroids.items():
        dist = np.linalg.norm(x_target - centroid)
        distances[pop] = float(dist)

    # Convert distances to probabilities (softmax-like)
    if distances:
        min_dist = min(distances.values())
        if min_dist < 1e-10:
            min_dist = 1e-10
        weights = {p: 1.0 / max(d, 1e-10) for p, d in distances.items()}
        total = sum(weights.values())
        probs = {p: round(w / total, 4) for p, w in weights.items()}
    else:
        probs = {p: 1.0 / len(pops) for p in pops}

    best = min(distances, key=distances.get)
    return {
        "assigned_population": best,
        "probabilities": probs,
        "centroid_distances": {p: round(d, 4) for p, d in distances.items()},
    }


def classify_knn(X_ref: np.ndarray, x_target: np.ndarray, y_ref: np.ndarray,
                 k: int = 15) -> Dict:
    """k-NN classification in PCA space (weighted by 1/distance)."""
    # Compute Euclidean distances to all reference samples
    diffs = X_ref - x_target
    distances = np.sqrt(np.sum(diffs ** 2, axis=1))
    # Get k nearest
    k = min(k, len(distances))
    top_k = np.argpartition(distances, k)[:k]
    top_distances = distances[top_k]
    top_indices = top_k

    votes = {}
    for i in range(len(top_indices)):
        d = top_distances[i]
        idx = top_indices[i]
        pop = str(y_ref[idx])
        weight = 1.0 / max(d, 1e-10)
        votes[pop] = votes.get(pop, 0.0) + weight

    total = sum(votes.values())
    probs = {p: round(v / total, 4) for p, v in sorted(votes.items(), key=lambda x: -x[1])}

    best = max(votes, key=votes.get)
    return {
        "assigned_population": best,
        "probabilities": probs,
    }


def ensemble(knn_result: Dict, centroid_result: Dict) -> Dict:
    """Combine k-NN and centroid results with equal weight."""
    combined = {}
    all_pops = set(list(knn_result["probabilities"].keys()) + list(centroid_result["probabilities"].keys()))
    for pop in all_pops:
        knn_p = knn_result["probabilities"].get(pop, 0.0)
        cent_p = centroid_result["probabilities"].get(pop, 0.0)
        combined[pop] = (knn_p + cent_p) / 2.0

    total = sum(combined.values())
    if total > 0:
        combined = {p: round(v / total, 4) for p, v in combined.items()}

    best = max(combined, key=combined.get)
    best_prob = combined[best]

    # Confidence based on max probability and method agreement
    knn_best = knn_result["assigned_population"]
    cent_best = centroid_result["assigned_population"]
    methods_agree = knn_best == cent_best

    if best_prob >= 0.90 and methods_agree:
        confidence = "HIGH"
    elif best_prob >= 0.70 and methods_agree:
        confidence = "MODERATE"
    elif best_prob >= 0.50:
        confidence = "LOW"
    else:
        confidence = "LOW"

    return {
        "assigned_population": best,
        "confidence": confidence,
        "max_probability": round(best_prob, 4),
        "methods_agree": methods_agree,
        "knn_best": knn_best,
        "centroid_best": cent_best,
        "posterior_probabilities": combined,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Sub-Continental Ancestry Classifier (EUR fine-structure)")
    parser.add_argument("--ref-pcs", default="pca/1000G_pcs.eigenvec")
    parser.add_argument("--target-pcs", default="pca/target_pcs.eigenvec")
    parser.add_argument("--pop-panel", default="reference/1000G_full/population_panel.txt")
    parser.add_argument("--output-dir", "-o", default="pca")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    logger.info("═══ Sub-Continental Ancestry Classifier ═══")

    # Load data
    X_ref, x_target, y_ref, pc_cols = load_data(args.ref_pcs, args.target_pcs, args.pop_panel)

    # Classify
    knn_result = classify_knn(X_ref, x_target, y_ref, k=15)
    centroid_result = classify_centroid(X_ref, x_target, y_ref, EUR_SUB_POPS)
    result = ensemble(knn_result, centroid_result)

    logger.info(f"  Assigned: {result['assigned_population']} ({SUB_POP_NAMES.get(result['assigned_population'], '?')})")
    logger.info(f"  Confidence: {result['confidence']} (max_prob={result['max_probability']:.3f})")
    logger.info(f"  Methods agree: {result['methods_agree']} (k-NN={knn_result['assigned_population']}, centroid={centroid_result['assigned_population']})")

    # Build output with full sub-population details
    sub_pops_available = []
    for pop in EUR_SUB_POPS:
        sub_pops_available.append({
            "code": pop,
            "name": SUB_POP_NAMES.get(pop, pop),
            "description": SUB_POP_DESCRIPTIONS.get(pop, ""),
            "n_reference_samples": int((y_ref == pop).sum()),
        })

    output = {
        "assigned_sub_population": result["assigned_population"],
        "sub_population_name": SUB_POP_NAMES.get(result["assigned_population"], result["assigned_population"]),
        "sub_population_description": SUB_POP_DESCRIPTIONS.get(result["assigned_population"], ""),
        "confidence": result["confidence"],
        "max_probability": result["max_probability"],
        "methods_agree": result["methods_agree"],
        "posterior_probabilities": result["posterior_probabilities"],
        "knn_assignment": knn_result["assigned_population"],
        "centroid_assignment": centroid_result["assigned_population"],
        "n_reference_samples": len(y_ref),
        "sub_populations_available": sub_pops_available,
        "method": "PCA ensemble (k-NN + centroid distance) on EUR-only 1000G subset",
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
    }

    out_path = Path(args.output_dir) / "subcontinental_assignment.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(output, fh, indent=2)
    logger.info(f"  ✅ Output: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
