#!/usr/bin/env python3
"""
Compute permutation-based empirical p-values for drug reversal scores
using the actual signature-level reversal-score pool as the null.

Saves output to: Reversal/C4A_drug_reversals_significance.tsv
"""

import os
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# --------- USER ADJUSTABLE ---------
N_PERMUTATIONS = 5000  # increase for more precision (slower)
RANDOM_SEED = 0
# ----------------------------------

np.random.seed(RANDOM_SEED)

# input filenames (preserve your paths)
PERT_LEVEL_PATH = "Reversal/connectivity_by_pert.tsv"
SID_LEVEL_PATH  = "Reversal/connectivity_by_sid.tsv"
COL_META = "lincs_col_meta_parsed.tsv"
PERT_INFO = "GSE92742_Broad_LINCS_pert_info.txt"

# output
OUTPATH = "Reversal/C4A_drug_reversals_significance.tsv"

print("[DEBUG] Loading files...")

# Prefer pert-level file if it contains signature-level rows (sig_id), else fallback to sid file
if os.path.exists(PERT_LEVEL_PATH):
    pert_scores = pd.read_csv(PERT_LEVEL_PATH, sep="\t", dtype=str, low_memory=False)
    if "sig_id" not in pert_scores.columns and os.path.exists(SID_LEVEL_PATH):
        print("[WARN] connectivity_by_pert.tsv lacks 'sig_id' column — falling back to connectivity_by_sid.tsv")
        pert_scores = pd.read_csv(SID_LEVEL_PATH, sep="\t", dtype=str, low_memory=False)
else:
    if os.path.exists(SID_LEVEL_PATH):
        pert_scores = pd.read_csv(SID_LEVEL_PATH, sep="\t", dtype=str, low_memory=False)
    else:
        print("[FATAL] Neither connectivity_by_pert.tsv nor connectivity_by_sid.tsv found. Exiting.")
        sys.exit(1)

# load meta (used only if needed)
if os.path.exists(COL_META):
    col_meta = pd.read_csv(COL_META, sep="\t", dtype=str, low_memory=False)
else:
    col_meta = None

if os.path.exists(PERT_INFO):
    pert_info = pd.read_csv(PERT_INFO, sep="\t", dtype=str, low_memory=False)
else:
    pert_info = None

print(f"[DEBUG] Loaded pert_scores with columns: {list(pert_scores.columns)}")

# Ensure necessary columns exist
required_cols = ["pert_id", "reversal_score"]
if "pert_id" not in pert_scores.columns:
    # Try to derive pert_id via col_meta if possible
    if col_meta is not None and "id" in col_meta.columns and "full_pert_batch" in col_meta.columns:
        print("[DEBUG] Deriving pert_id from full_pert_batch via sig_id mapping...")
        pert_scores = pert_scores.merge(col_meta[["id", "full_pert_batch"]].rename(columns={"id": "sig_id", "full_pert_batch": "pert_id"}),
                                       on="sig_id", how="left", suffixes=("", "_meta"))
        if "pert_id_meta" in pert_scores.columns:
            pert_scores["pert_id"] = pert_scores["pert_id"].fillna(pert_scores["pert_id_meta"])
            pert_scores.drop(columns=["pert_id_meta"], inplace=True)
    if "pert_id" not in pert_scores.columns:
        print("[FATAL] pert_id not found and could not derive. Please provide a signature-level file with 'pert_id'.")
        sys.exit(1)

if "reversal_score" not in pert_scores.columns:
    print("[FATAL] 'reversal_score' column not found in signature-level connectivity file.")
    sys.exit(1)

# force numeric
pert_scores["reversal_score"] = pd.to_numeric(pert_scores["reversal_score"], errors="coerce")
# drop rows with missing reversal_score or pert_id
pert_scores = pert_scores.dropna(subset=["pert_id", "reversal_score"]).copy()

print(f"[DEBUG] Signature-level rows available for null pool: {len(pert_scores)}")

# If pert_info provided, filter to trt_cp so null uses only compound signatures (optional)
if pert_info is not None and "pert_type" in pert_info.columns:
    cp_set = set(pert_info.loc[pert_info["pert_type"] == "trt_cp", "pert_id"].astype(str).values)
    if len(cp_set) > 0:
        before = len(pert_scores)
        pert_scores = pert_scores[pert_scores["pert_id"].astype(str).isin(cp_set)].copy()
        print(f"[DEBUG] Filtered signature pool to compound perturbations: {before} -> {len(pert_scores)}")

# Build observed aggregated table (per-pert)
agg_cols = []
if "spearman" in pert_scores.columns:
    agg_cols.append("spearman")
# observed mean reversal_score per pert
obs = (
    pert_scores.groupby(["pert_id"])
    .agg(
        reversal_score_obs = ("reversal_score", "mean"),
        n_signatures = ("reversal_score", "size")
    )
    .reset_index()
)

# If pert_iname exists in pert_scores, keep; else try to merge with pert_info
if "pert_iname" in pert_scores.columns:
    # take first pert_iname per pert_id
    names = pert_scores.groupby("pert_id").agg(pert_iname=("pert_iname", "first")).reset_index()
    obs = obs.merge(names, on="pert_id", how="left")
elif pert_info is not None and "pert_iname" in pert_info.columns:
    obs = obs.merge(pert_info[["pert_id", "pert_iname"]].drop_duplicates(), on="pert_id", how="left")
else:
    obs["pert_iname"] = obs["pert_id"]

print(f"[DEBUG] Number of unique perturbagens to test: {len(obs)}")

# Build the null pool: array of signature-level reversal scores
null_pool = pert_scores["reversal_score"].values
pool_size = len(null_pool)
if pool_size == 0:
    print("[FATAL] No signature-level reversal scores available to build null distribution.")
    sys.exit(1)

# Pre-allocate arrays for results
emp_pvals = np.zeros(len(obs), dtype=float)
z_scores  = np.zeros(len(obs), dtype=float)

print(f"[INFO] Running {N_PERMUTATIONS} permutations per drug (this can take a while)...")

# To speed: create a large (N_permutations x max_k) random index matrix for sampling (with replacement),
# but we'll sample per drug because k varies.
for i, row in tqdm(obs.iterrows(), total=len(obs), desc="per-drug perms"):
    k = int(row["n_signatures"])
    observed_mean = float(row["reversal_score_obs"])
    # sampling strategy: if pool_size >= k and k>1, sample without replacement; else with replacement
    if pool_size >= k:
        # draw N_PERMUTATIONS samples of k indices without replacement
        # For efficiency, sample indices for each perm separately
        perm_means = []
        for _ in range(N_PERMUTATIONS):
            idx = np.random.choice(pool_size, size=k, replace=False)
            perm_means.append(null_pool[idx].mean())
        perm_means = np.array(perm_means)
    else:
        # pool smaller than k -> sample with replacement
        perm_means = np.random.choice(null_pool, size=(N_PERMUTATIONS, k), replace=True).mean(axis=1)

    # empirical one-sided p-value: proportion of permuted means >= observed_mean
    # add +1 numerator and denominator for conservative estimate
    emp_p = (np.sum(perm_means >= observed_mean) + 1) / (N_PERMUTATIONS + 1)
    emp_pvals[i] = emp_p

    # z-score relative to permuted means
    mu = perm_means.mean()
    sigma = perm_means.std(ddof=1) if perm_means.std(ddof=1) > 0 else np.nan
    if np.isfinite(sigma) and sigma > 0:
        z_scores[i] = (observed_mean - mu) / sigma
    else:
        z_scores[i] = np.nan

# Add results to obs
obs["empirical_pval"] = emp_pvals
obs["z_score"] = z_scores

# Multiple testing correction (Benjamini-Hochberg)
obs["FDR"] = multipletests(obs["empirical_pval"], method="fdr_bh")[1]
obs["significant"] = obs["FDR"] < 0.05

# Order by reversal score
obs = obs.sort_values("reversal_score_obs", ascending=False).reset_index(drop=True)

# Reformat columns for output
out = obs.rename(columns={"reversal_score_obs": "reversal_score", "n_signatures": "n_signatures"})
cols_out = ["pert_id", "pert_iname", "reversal_score", "z_score", "n_signatures", "empirical_pval", "FDR", "significant"]
# ensure all exist
for c in cols_out:
    if c not in out.columns:
        out[c] = np.nan

out = out[cols_out]

# Save
out.to_csv(OUTPATH, sep="\t", index=False)
print(f"[SUCCESS] Wrote: {OUTPATH}")
print(f"[INFO] Significant hits (FDR < 0.05): {out['significant'].sum()} / {len(out)}")

# Quick top-20 print
print("\n=== Top 20 (by reversal_score) with significance ===")
print(out.head(20).to_string(index=False))
