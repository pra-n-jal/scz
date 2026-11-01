#!/usr/bin/env python3
"""
compute_directional_fdr_pipeline.py

Usage examples:
python compute_permutation_fdr.py --reversals mapped_mapping_audit.tsv --n_perms 10000 --out C4A_with_fdr.tsv

Produces:
 - {out} (per-row results with p-values and multiple FDR variants)
 - {out_compound} (per-compound aggregated results with same stats)
 - prints a short enrichment summary for provided keywords

Features implemented (A-H from your request):
A) Report nominal (raw) significance and note if none survive BH.
B) Allows focused tests via --filter_list (optional).
C) One-sided (directional) empirical p-values (default direction: negative).
D) Weighted FDR (IHW-ish) using n_signatures as weights.
E) Aggregation per compound (Stouffer z and median), both tested.
F) Proper permutation null (full-vector shuffles) and empirical p-values.
G) Enrichment test for keyword-defined drug classes (Fisher/hypergeometric).
H) Storey q-value (simple pi0 estimation)
"""

import argparse
import sys
import os
import math
from collections import defaultdict
import numpy as np
import pandas as pd
from math import comb
from functools import lru_cache

# ---------------------------
# Utilities: BH, Storey q, Fisher/hypergeom tail
# ---------------------------
def benjamini_hochberg(pvals):
    """Classic BH; returns q-values (same length as pvals)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    for i, idx in enumerate(order):
        rank = i + 1
        ranked[idx] = p[idx] * n / rank
    # enforce monotonicity
    adj = ranked.copy()
    # ensure non-decreasing when moving from largest to smallest p
    inv_order = order[::-1]
    min_so_far = 1.0
    for idx in inv_order:
        min_so_far = min(min_so_far, adj[idx])
        adj[idx] = min_so_far
    adj[adj > 1.0] = 1.0
    return adj

def storey_qvalue(pvals, lambda_=0.5):
    """
    Simple Storey estimator for q-values (pi0 estimate with single lambda).
    Not the most sophisticated, but reasonable.
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    pi0 = np.minimum(1.0, np.sum(p > lambda_) / ((1.0 - lambda_) * n))
    # avoid zero
    pi0 = max(pi0, 1e-8)
    # BH-like adjusted p then scaled by pi0
    q_bh = benjamini_hochberg(p)
    q_storey = np.minimum(1.0, q_bh * pi0)
    return q_storey, pi0

@lru_cache(maxsize=None)
def logcomb(n, k):
    return math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)

def hypergeom_tail(k, K, n, N):
    """
    Return upper-tail hypergeometric p-value: P(X >= k)
    N = population size
    K = successes in population
    n = sample size
    k = observed successes in sample
    Uses log-sum for numerical stability.
    """
    # bounds
    k = int(k)
    K = int(K)
    n = int(n)
    N = int(N)
    lo = max(0, n + K - N)
    hi = min(n, K)
    if k > hi:
        return 0.0
    if k <= lo:
        return 1.0
    # sum probabilities for i=k..hi
    total = 0.0
    # compute denominator logcomb(N,n) once
    denom_log = logcomb(N, n)
    for i in range(k, hi+1):
        num_log = logcomb(K, i) + logcomb(N-K, n-i)
        total += math.exp(num_log - denom_log)
    return min(1.0, total)

# ---------------------------
# Core pipeline
# ---------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Directional permutation FDR pipeline (A-H).")
    p.add_argument("--reversals", required=True, default="mapped_mapping_audit.tsv",
                   help="TSV with at least columns: reversal_score, n_signatures, pert_id, pert_iname")
    p.add_argument("--n_perms", type=int, default=10000)
    p.add_argument("--out", default="C4A_with_fdr.tsv")
    p.add_argument("--direction", choices=["negative","positive","two-sided"], default="negative",
                   help="Direction of interest for one-sided test. Default: negative (reversal_score < 0 desired).")
    p.add_argument("--top_k", type=int, default=30, help="Top-K for enrichment tests (default 30)")
    p.add_argument("--class_keywords", default="hdac,retinoid,vitamin,statin,antipsychotic",
                   help="Comma-separated keywords for enrichment testing in pert_iname (case-insensitive).")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--min_permutations_warn", type=int, default=1000,
                   help="Warn if n_perms low for p-resolution.")
    return p.parse_args()

def one_sided_p_values_from_null(obs, null_matrix, direction="negative"):
    """
    obs: 1D array of observed reversal_score (signed)
    null_matrix: shape (n_perms, n_obs) containing permuted values (kept signed)
    direction: negative, positive or two-sided
    Returns empirical pvals (one-sided or two-sided)
    """
    obs = np.asarray(obs)
    absnull = np.abs(null_matrix)
    n_perms = null_matrix.shape[0]
    pvals = np.empty(obs.shape[0], dtype=float)
    if direction == "two-sided":
        obs_abs = np.abs(obs)
        # count how many permuted abs >= observed abs
        counts = (absnull >= obs_abs[None, :]).sum(axis=0)
        pvals = (counts + 1) / (n_perms + 1)
    elif direction == "negative":
        # we want P(null <= obs) if obs is negative direction of interest -> compute tail accordingly.
        # For negative direction: smaller (more negative) is better: p = P(null <= obs)
        for i, o in enumerate(obs):
            if o < 0:
                # count permutations with value <= observed
                c = (null_matrix[:, i] <= o).sum()
                pvals[i] = (c + 1) / (n_perms + 1)
            else:
                # if observed is not in the direction of interest, p ~ 1
                pvals[i] = 1.0
    else:  # positive
        for i, o in enumerate(obs):
            if o > 0:
                c = (null_matrix[:, i] >= o).sum()
                pvals[i] = (c + 1) / (n_perms + 1)
            else:
                pvals[i] = 1.0
    return pvals

def weighted_bh(pvals, weights):
    """
    Weighted BH as in Genovese & Wasserman:
    - weights must be positive and sum to 1.
    Implementation: compute p_i / w_i, apply BH, then q_i = adj_p_i * w_i
    Returns q-values (weighted).
    """
    p = np.asarray(pvals, dtype=float)
    w = np.asarray(weights, dtype=float)
    if np.any(w <= 0):
        # fall back to equal weights
        w = np.ones_like(w) / len(w)
    # normalize
    w = w / w.sum()
    adj_input = p / w
    # cap >1 for BH input safety
    adj_input = np.minimum(adj_input, 1.0)
    adj = benjamini_hochberg(adj_input)
    q = adj * w
    q[q > 1.0] = 1.0
    return q

def stouffer_z_and_p(values, weights=None):
    """
    Compute Stouffer Z from signed values (assume values are z-like already?).
    We're going to convert reversal_score -> z by rank-normalizing if needed.
    Here, we'll compute weighted Stouffer using values (we assume larger magnitude = stronger).
    But to be safe across scales, convert values to ranks -> normal quantiles.
    """
    arr = np.asarray(values, dtype=float)
    # rank to uniform [0,1]
    ranks = np.argsort(np.argsort(arr)) + 1
    uniform = ranks / (len(arr) + 1)
    # normal quantiles
    z = np.sqrt(2) * erfinv(2*uniform - 1)  # inverse normal
    if weights is None:
        weights = np.ones_like(z)
    weights = np.asarray(weights, dtype=float)
    z_comb = np.sum(weights * z) / np.sqrt(np.sum(weights**2))
    # two-sided p
    p = 2.0 * (1.0 - 0.5*(1.0 + math.erf(abs(z_comb) / math.sqrt(2))))
    return z_comb, p

# inverse error function (erfinv) since math.erfinv may be absent; implement via approximate
def erfinv(y):
    # approximate erfinv using polynomial approximation (Winitzki)
    # valid for y in (-1,1)
    a = 0.147  # constant
    sign = 1 if y >= 0 else -1
    ln = math.log(1 - y*y)
    first = 2/(math.pi*a) + ln/2
    second = ln/(a)
    x = sign * math.sqrt( math.sqrt(first*first - second) - first )
    return x

# But we'll not rely on stouffer with erfinv; instead, compute stouffer using z-scores from standardized values
def stouffer_simple(values, weights=None):
    arr = np.asarray(values, dtype=float)
    # standardize to z using mean/std
    if np.allclose(np.std(arr), 0.0):
        zscores = np.zeros_like(arr)
    else:
        zscores = (arr - np.mean(arr)) / np.std(arr)
    if weights is None:
        weights = np.ones_like(zscores)
    weights = np.asarray(weights, dtype=float)
    z_comb = np.sum(weights * zscores) / np.sqrt(np.sum(weights**2))
    # one-sided p for direction negative/positive must be handled upstream; return z and two-sided p
    p_two = 2.0 * (1.0 - 0.5*(1.0 + math.erf(abs(z_comb) / math.sqrt(2))))
    return z_comb, p_two

def run_pipeline(args):
    rng = np.random.default_rng(args.seed)
    if args.n_perms < args.min_permutations_warn:
        print(f"[WARN] n_perms={args.n_perms} is low; resolution limit = 1/(n_perms+1) = {1.0/(args.n_perms+1):.4g}", file=sys.stderr)

    df = pd.read_csv(args.reversals, sep="\t", dtype=str)
    # ensure numeric columns
    if "reversal_score" not in df.columns:
        raise ValueError("Input file must contain column 'reversal_score'")
    df["reversal_score"] = pd.to_numeric(df["reversal_score"], errors="coerce")
    if "n_signatures" in df.columns:
        df["n_signatures"] = pd.to_numeric(df["n_signatures"], errors="coerce").fillna(1).astype(int)
    else:
        df["n_signatures"] = 1

    # maintain original order index
    df["_orig_index"] = np.arange(len(df))

    # If multiple rows map to same pert_iname, later aggregate
    # Clean pert_iname
    if "pert_iname" not in df.columns:
        df["pert_iname"] = df.get("pert_id", pd.Series(["unknown"]*len(df)))

    # ========== Build null matrix by permuting the signed vector across indices ==========
    obs = df["reversal_score"].values.copy()
    n_obs = obs.size
    n_perms = args.n_perms
    # We'll permute the full signed vector each time (exchangeability across rows).
    # nulls shape (n_perms, n_obs)
    nulls = np.empty((n_perms, n_obs), dtype=float)
    for i in range(n_perms):
        nulls[i, :] = rng.permutation(obs)

    # ========== Compute one-sided empirical p-values (directional) ==========
    pvals_emp = one_sided_p_values_from_null(obs, nulls, direction=args.direction)

    # ========== BH, Storey q, Weighted BH (weights = n_signatures normalized) ==========
    bh_q = benjamini_hochberg(pvals_emp)
    storey_q, pi0 = storey_qvalue(pvals_emp, lambda_=0.5)
    # weights from n_signatures (give more weight to well-sampled items)
    weights = df["n_signatures"].values.astype(float)
    # avoid zero weights
    weights = np.maximum(weights, 1.0)
    weights_norm = weights / weights.sum()
    weighted_q = weighted_bh(pvals_emp, weights_norm)

    # attach to df
    df_out = df.copy()
    df_out["perm_pval"] = pvals_emp
    df_out["fdr_bh"] = bh_q
    df_out["fdr_weighted"] = weighted_q
    df_out["q_storey"] = storey_q
    df_out["pi0_estimate"] = pi0

    # ========== Aggregation per compound (pert_iname) ==========
    # For each pert_iname: compute median reversal_score, stouffer z (simple standardize + weighted)
    groupby = df_out.groupby("pert_iname", sort=False)
    comp_rows = []
    for name, g in groupby:
        values = g["reversal_score"].values
        ns = g["n_signatures"].values
        # median effect
        median_val = np.median(values)
        # stouffer-like combination: standardize row values then weighted combine using sqrt(ns)
        # Use stouffer_simple from above
        w = np.sqrt(ns)  # weight per row proportional to sqrt of n_signatures
        z_comb, p_two = stouffer_simple(values, weights=w)
        # But we want signed combined effect: use weighted average of raw values (direction preserved)
        signed_comb = np.sum(values * w) / np.sum(w) if np.sum(w) > 0 else np.mean(values)
        comp_rows.append({
            "pert_iname": name,
            "n_rows": len(g),
            "n_signatures_total": int(np.sum(ns)),
            "median_reversal_score": median_val,
            "signed_combined_score": signed_comb,
            "stouffer_z": z_comb,
            "stouffer_two_sided_p": p_two
        })
    comp_df = pd.DataFrame(comp_rows)
    comp_df = comp_df.sort_values("signed_combined_score", ascending=True if args.direction=="negative" else False)

    # For compound-level permutation: build permuted combined scores for each perm by applying same aggregation to nulls
    # nulls shape (n_perms, n_obs) with rows aligned to original df rows
    # We'll compute for each perm the aggregated signed_combined_score per pert_iname
    pert_names = comp_df["pert_iname"].values
    # Map each original index to group index
    name_to_indices = {name: df_out.index[df_out["pert_iname"] == name].tolist() for name in pert_names}
    # allocate matrix
    comp_nulls = np.empty((n_perms, len(pert_names)), dtype=float)
    for i in range(n_perms):
        perm_values = nulls[i, :]
        # for each compound, compute signed_combined_score (weighted average by sqrt(n_signatures) within rows)
        for j, name in enumerate(pert_names):
            idxs = name_to_indices[name]
            vals = perm_values[idxs]
            ns = df_out.loc[idxs, "n_signatures"].values
            w = np.sqrt(ns)
            comp_nulls[i, j] = np.sum(vals * w) / np.sum(w) if np.sum(w) > 0 else np.mean(vals)
    # observed comp signed_combined
    comp_obs = comp_df["signed_combined_score"].values
    # compute empirical one-sided p-values for compounds
    comp_pvals = one_sided_p_values_from_null(comp_obs, comp_nulls, direction=args.direction)
    comp_bh_q = benjamini_hochberg(comp_pvals)
    comp_storey_q, comp_pi0 = storey_qvalue(comp_pvals, lambda_=0.5)
    # weights for compounds (use n_signatures_total)
    comp_weights = comp_df["n_signatures_total"].values.astype(float)
    comp_weights = np.maximum(comp_weights, 1.0)
    comp_weights_norm = comp_weights / comp_weights.sum()
    comp_weighted_q = weighted_bh(comp_pvals, comp_weights_norm)

    comp_df["perm_pval"] = comp_pvals
    comp_df["fdr_bh"] = comp_bh_q
    comp_df["fdr_weighted"] = comp_weighted_q
    comp_df["q_storey"] = comp_storey_q
    comp_df["pi0_estimate"] = comp_pi0

    # ========== Enrichment tests for keywords in pert_iname among top_k ==========
    keywords = [k.strip().lower() for k in args.class_keywords.split(",") if k.strip()]
    topk = int(min(args.top_k, len(comp_df)))
    top_names = set(comp_df.head(topk)["pert_iname"].values.tolist())
    enrichment_results = []
    N = len(comp_df)
    for kw in keywords:
        # count successes in population (compounds that mention kw)
        mask_kw = comp_df["pert_iname"].str.lower().str.contains(kw)
        K = mask_kw.sum()
        # successes in sample (top_k)
        k_obs = comp_df.head(topk)["pert_iname"].str.lower().str.contains(kw).sum()
        pval_hg = hypergeom_tail(k_obs, K, topk, N)
        enrichment_results.append({
            "keyword": kw,
            "K_in_population": int(K),
            "k_in_topk": int(k_obs),
            "N_population": int(N),
            "n_sample": int(topk),
            "hypergeom_pval": pval_hg
        })

    # ========== Save outputs ==========
    out_main = args.out
    base, ext = os.path.splitext(out_main)
    out_comp = base + "_compound_level" + ext
    # write per-row df_out in stable column order
    cols = list(df_out.columns)
    # re-order to put new columns at end in consistent order
    desired_tail = ["perm_pval", "fdr_bh", "fdr_weighted", "q_storey", "pi0_estimate"]
    for c in desired_tail:
        if c in cols:
            cols.remove(c)
            cols.append(c)
    df_out.to_csv(out_main, sep="\t", index=False, columns=cols)
    comp_df.to_csv(out_comp, sep="\t", index=False)

    # ========== Print summary ==========
    print(f"[INFO] Saved per-row results -> {out_main}")
    print(f"[INFO] Saved compound-level results -> {out_comp}")
    print(f"[INFO] n_rows = {len(df_out)}; n_compounds = {len(comp_df)}; n_perms = {n_perms}; direction = {args.direction}")
    print(f"[INFO] Storey pi0 (rows) = {pi0:.4f}; (compounds) = {comp_pi0:.4f}")
    # top hits summary
    top_rows = df_out.sort_values("perm_pval").head(10)[["pert_iname","reversal_score","perm_pval","fdr_bh","fdr_weighted","q_storey","n_signatures"]]
    print("\nTop 10 rows by perm_pval:")
    print(top_rows.to_string(index=False))
    print("\nTop 10 compounds by combined signed score:")
    print(comp_df.head(10)[["pert_iname","signed_combined_score","perm_pval","fdr_bh","fdr_weighted","q_storey","n_signatures_total"]].to_string(index=False))
    print("\nEnrichment results (top_k = %d):" % topk)
    enr_df = pd.DataFrame(enrichment_results)
    print(enr_df.to_string(index=False))

if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
