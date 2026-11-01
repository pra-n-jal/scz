#!/usr/bin/env python3
"""
Compute permutation-based empirical p-values for drug reversal scores
using the actual signature-level reversal-score pool as the null.
"""

import os
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Compute permutation-based empirical p-values for drug reversal scores.")
    parser.add_argument("--pert-scores", default="results/connectivity_by_pert.tsv",
                        help="Path to connectivity_by_pert.tsv")
    parser.add_argument("--sid-scores", default="results/connectivity_by_sid.tsv",
                        help="Path to connectivity_by_sid.tsv")
    parser.add_argument("--col-meta", default="lincs_col_meta_parsed.tsv",
                        help="Path to parsed LINCS column metadata")
    parser.add_argument("--pert-info", default="data/raw/GSE92742_Broad_LINCS_pert_info.txt",
                        help="Path to GSE92742_Broad_LINCS_pert_info.txt")
    parser.add_argument("--out", default="results/C4A_drug_reversals_with_significance.tsv",
                        help="Output TSV path")
    parser.add_argument("--n-permutations", type=int, default=5000,
                        help="Number of permutations (default: 5000)")
    parser.add_argument("--seed", type=int, default=0,
                        help="Random seed (default: 0)")
    return parser.parse_args()


def main():
    args = parse_args()
    np.random.seed(args.seed)

    pert_path = args.pert_scores
    if not os.path.exists(pert_path):
        alt = "Reversal/connectivity_by_pert.tsv"
        if os.path.exists(alt):
            pert_path = alt

    sid_path = args.sid_scores
    if not os.path.exists(sid_path):
        alt = "Reversal/connectivity_by_sid.tsv"
        if os.path.exists(alt):
            sid_path = alt

    print("[DEBUG] Loading files...")
    if os.path.exists(pert_path):
        pert_scores = pd.read_csv(pert_path, sep="\t", dtype=str, low_memory=False)
        if "sig_id" not in pert_scores.columns and os.path.exists(sid_path):
            print("[WARN] connectivity_by_pert.tsv lacks 'sig_id' column — falling back to connectivity_by_sid.tsv")
            pert_scores = pd.read_csv(sid_path, sep="\t", dtype=str, low_memory=False)
    else:
        if os.path.exists(sid_path):
            pert_scores = pd.read_csv(sid_path, sep="\t", dtype=str, low_memory=False)
        else:
            print("[FATAL] Neither connectivity_by_pert.tsv nor connectivity_by_sid.tsv found. Exiting.")
            sys.exit(1)

    col_meta = pd.read_csv(args.col_meta, sep="\t", dtype=str, low_memory=False) if os.path.exists(args.col_meta) else None
    pert_info = pd.read_csv(args.pert_info, sep="\t", dtype=str, low_memory=False) if os.path.exists(args.pert_info) else None

    print(f"[DEBUG] Loaded pert_scores with columns: {list(pert_scores.columns)}")

    # Score column detection
    score_col = None
    for cand in ["reversal_score", "spearman"]:
        if cand in pert_scores.columns:
            score_col = cand
            break
    if score_col is None:
        print("[FATAL] Neither 'reversal_score' nor 'spearman' column present.")
        sys.exit(1)

    pert_scores[score_col] = pd.to_numeric(pert_scores[score_col], errors="coerce")
    valid_scores = pert_scores[score_col].dropna().values
    if len(valid_scores) == 0:
        print("[FATAL] No numeric scores available to build null distribution.")
        sys.exit(1)

    # Observed aggregation by pert_iname or pert_id
    id_col = "pert_iname" if "pert_iname" in pert_scores.columns else ("pert_id" if "pert_id" in pert_scores.columns else pert_scores.columns[0])
    obs = pert_scores.groupby(id_col)[score_col].agg(["mean", "count"]).reset_index()
    obs = obs.rename(columns={"mean": "reversal_score_obs", "count": "n_signatures"})
    obs = obs[obs["n_signatures"] >= 1].copy()

    null_pool = valid_scores
    null_mean = np.mean(null_pool)
    null_std = np.std(null_pool)

    # Permutation test per compound
    p_vals = []
    z_scores = []
    print(f"[INFO] Computing empirical p-values across {len(obs)} compounds ({args.n_permutations} permutations)...")

    for _, row in tqdm(obs.iterrows(), total=len(obs)):
        k = int(row["n_signatures"])
        obs_val = row["reversal_score_obs"]
        if k == 1:
            emp_p = (np.sum(null_pool >= obs_val) + 1.0) / (len(null_pool) + 1.0)
            z = (obs_val - null_mean) / (null_std if null_std > 0 else 1.0)
        else:
            perm_means = np.random.choice(null_pool, size=(args.n_permutations, k), replace=True).mean(axis=1)
            emp_p = (np.sum(perm_means >= obs_val) + 1.0) / (args.n_permutations + 1.0)
            perm_std = np.std(perm_means)
            z = (obs_val - np.mean(perm_means)) / (perm_std if perm_std > 0 else 1.0)

        p_vals.append(emp_p)
        z_scores.append(z)

    obs["empirical_pval"] = p_vals
    obs["z_score"] = z_scores

    # FDR correction
    reject, q_vals, _, _ = multipletests(obs["empirical_pval"], method="fdr_bh")
    obs["FDR"] = q_vals
    obs["significant"] = reject
    obs = obs.sort_values("reversal_score_obs", ascending=False).reset_index(drop=True)
    obs = obs.rename(columns={"reversal_score_obs": "reversal_score"})

    out_p = Path(args.out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    obs.to_csv(out_p, sep="\t", index=False)
    print(f"[SUCCESS] Wrote: {out_p}")
    print(f"[INFO] Significant hits (FDR < 0.05): {obs['significant'].sum()} / {len(obs)}")


if __name__ == "__main__":
    main()
