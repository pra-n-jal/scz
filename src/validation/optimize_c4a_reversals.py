#!/usr/bin/env python3
"""
optimize_c4a_reversals.py

Usage:
    python optimize_c4a_reversals.py \
        --rev C4A_drug_reversals.tsv \
        --pert pert_info.csv \
        --out_prefix C4A_drug_reversals.optimized \
        [--permute drug_gene_matrix.tsv c4a_signature.tsv] \
        [--topN_permute 200] [--n_perm 1000]

Outputs (prefix-based):
 - {prefix}.cleaned.tsv                 cleaned input (duplicate header rows removed)
 - {prefix}.merged.tsv                  merged with pert_info columns
 - {prefix}.unmapped_pertids.tsv        unmapped pert_ids for manual curation (prioritized)
 - {prefix}.aggregated_by_name.tsv      aggregated metrics per pert_iname
 - {prefix}.bias_adjusted.tsv           contains adjusted scores (log-n downweight, resid)
 - {prefix}.final_prioritized.tsv       final ranking with CNS priority and final_score
 - optional: {prefix}.permutation_top{N}.tsv  (if signature matrices provided)

Notes:
 - The script assumes the reversal file is tab-separated (TSV) with header:
   pert_id\tpert_iname\tspearman\treversal_score\tsingle_gene_value\tn_signatures
 - pert_info.csv must have header like:
   pert_id,pert_iname,pert_type,is_touchstone,inchi_key_prefix,inchi_key,canonical_smiles,pubchem_cid
 - RDKit: used only if canonical_smiles present and RDKit is installed; optional.
"""

import argparse, sys, os, re
import pandas as pd
import numpy as np

# Optional imports
try:
    import statsmodels.api as sm
except Exception:
    sm = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False

ONCOLOGY_KEYWORDS = [
    'cancer','carcinoma','tumor','onc','-tinib','chemotherapy','antineoplastic',
    'carboplatin','cisplatin','doxorubicin','cyclophosphamide','gemcitabine'
]

def is_probably_oncology(name):
    if pd.isna(name): return False
    s = name.lower()
    return any(k in s for k in ONCOLOGY_KEYWORDS)

def clean_input_rev(path):
    # read robustly (detect sep)
    df = pd.read_csv(path, sep=None, engine='python', dtype=str)
    # Remove accidental duplicate header rows where pert_id == 'pert_id'
    df = df[df['pert_id'].str.lower() != 'pert_id'].copy()
    # Strip whitespace
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # numeric coercion
    for col in ['reversal_score','single_gene_value','n_signatures','spearman']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def load_pert_info(path):
    # CSV with comma sep
    p = pd.read_csv(path, sep=',', dtype=str, low_memory=False)
    # standardize column names
    p.columns = [c.strip() for c in p.columns]
    return p

def flag_unmapped(row):
    name = row.get('pert_iname','') or ''
    pid = row.get('pert_id','') or ''
    # if pert_iname missing or equals pert_id -> unmapped
    if pd.isna(name) or name.strip()=='':
        return True
    if name.strip().upper() == pid.strip().upper():
        return True
    if re.match(r'^(BRD|CMAP|CPC|CRCGN)-', name.strip(), flags=re.I):
        return True
    # if no canonical_smiles or pubchem id in merged step -> will be flagged later
    return False

def aggregate_by_name(df):
    gp = df.groupby('pert_iname', dropna=False)
    agg = gp.agg(
        n_pert_ids = ('pert_id', lambda x: x.nunique()),
        total_signatures = ('n_signatures', 'sum'),
        median_score = ('reversal_score', 'median'),
        mean_score = ('reversal_score', 'mean'),
        n_rows = ('pert_id','count'),
        any_unmapped = ('unmapped_flag','any')
    ).reset_index()
    # Fill NaN numeric with 0
    agg[['total_signatures','n_rows']] = agg[['total_signatures','n_rows']].fillna(0).astype(int)
    return agg

def adjust_for_library_size(df_agg, method='log'):
    # method options: log, power
    df = df_agg.copy()
    df['log_n'] = np.log1p(df['total_signatures'])
    if method == 'log':
        df['adj_score'] = df['median_score'] / df['log_n'].replace(0, 1.0)
    elif method == 'power':
        df['adj_score'] = df['median_score'] / (np.power(df['total_signatures'].replace(0,1), 0.25))
    else:
        raise ValueError("method must be 'log' or 'power'")
    return df

def residualize_scores(df_agg):
    if sm is None:
        print("[WARN] statsmodels not available; skipping residualization.")
        df_agg['score_resid'] = np.nan
        return df_agg
    df = df_agg.copy()
    # create oncology flag based on pert_iname text heuristic
    df['onc_flag'] = df['pert_iname'].apply(lambda x: 1 if is_probably_oncology(x) else 0)
    X = pd.DataFrame({
        'log_n': np.log1p(df['total_signatures']),
        'onc': df['onc_flag']
    })
    X = sm.add_constant(X)
    y = df['median_score'].fillna(0.0)
    model = sm.OLS(y, X).fit()
    df['score_resid'] = model.resid
    df['resid_model_summary'] = model.rsquared
    return df

def bbb_heuristic_from_smiles(smiles):
    # RDKit-based simple scoring 0-1. 1 = very BBB-feasible (heuristic)
    if not RDKit_AVAILABLE:
        return np.nan
    if not isinstance(smiles, str) or smiles.strip()=='':
        return 0.0
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return 0.0
        mw = Descriptors.MolWt(m)
        psa = Descriptors.TPSA(m)
        hbd = Descriptors.NumHDonors(m)
        clogp = Descriptors.MolLogP(m)
    except Exception:
        return 0.0
    score = 0
    if mw < 450: score += 1
    if psa < 90: score += 1
    if hbd <= 3: score += 1
    if -1 <= clogp <= 4: score += 1
    return score / 4.0

def compute_cns_priority(df, require_smiles=False):
    d = df.copy()
    # approved_flag maybe present in pert_info; try to detect
    d['approved_flag'] = d['approved_flag'] if 'approved_flag' in d.columns else 0.0
    d['approved_flag'] = d['approved_flag'].fillna(0).astype(float)

    # If canonical_smiles exists, compute bbb_prob
    if 'canonical_smiles' in d.columns:
        d['bbb_prob'] = d['canonical_smiles'].fillna('').apply(bbb_heuristic_from_smiles)
        # If RDKit not available, bbb_prob will be NaN -> set to 0
        if not RDKit_AVAILABLE:
            d['bbb_prob'] = 0.0
    else:
        d['bbb_prob'] = 0.0
    # cytotoxic heuristic
    d['cytotoxic_flag'] = d['pert_iname'].fillna('').str.lower().apply(
        lambda s: 1 if any(k in s for k in ['cytotoxic','antineoplastic','chemotherapy','carboplatin','cisplatin','doxorubicin']) else 0
    )
    # Compose priority
    d['CNS_priority'] = 0.4 * d['approved_flag'] + 0.3 * d['bbb_prob'] + 0.2 * (1 - d['cytotoxic_flag'])
    # final_score uses adj_score if present, else median_score
    base = 'adj_score' if 'adj_score' in d.columns else 'median_score'
    d['final_score'] = d[base].fillna(0) * (1.0 + d['CNS_priority'].fillna(0))
    return d

def do_permutation_test(drug_mat_path, c4a_sig_path, topN=200, n_perm=1000, out_prefix='perm'):
    """
    Expects:
     - drug_gene_matrix.tsv: rows = signature ids (or pert_ids), columns = gene symbols, values = LFC or z
     - c4a_signature.tsv: two-column TSV: gene <tab> value (index = gene)
    This function will compute Spearman obs for each signature, then permute c4a vector to compute empirical p-values.
    Writes topN permutation results to file.
    """
    import time
    from scipy.stats import spearmanr
    print("[PERM] Loading matrices...")
    drug_mat = pd.read_csv(drug_mat_path, sep='\t', index_col=0)
    c4a = pd.read_csv(c4a_sig_path, sep='\t', index_col=0, header=None).iloc[:,0]
    common = drug_mat.columns.intersection(c4a.index)
    if len(common) < 50:
        print("[PERM] WARNING: small overlap of genes:", len(common))
    drug_mat = drug_mat[common]
    c4a = c4a[common]
    # compute observed
    obs = drug_mat.apply(lambda row: spearmanr(row.values, c4a.values).correlation, axis=1).rename('obs_r')
    # pick top absolute hits
    top_idx = obs.abs().sort_values(ascending=False).head(topN).index
    obs_top = obs.loc[top_idx]
    null_dict = {i: [] for i in top_idx}
    rng = np.random.default_rng(seed=12345)
    start = time.time()
    for i in range(n_perm):
        if (i+1) % 50 == 0:
            print(f"[PERM] iter {i+1}/{n_perm}")
        shuffled = rng.permutation(c4a.values)
        for sig in top_idx:
            r = spearmanr(drug_mat.loc[sig].values, shuffled).correlation
            null_dict[sig].append(r)
    # compute empirical p-values
    rows = []
    for sig in top_idx:
        obs_r = obs_top.loc[sig]
        null_arr = np.array(null_dict[sig])
        p_emp = (np.sum(null_arr <= obs_r) + 1) / (len(null_arr) + 1)
        rows.append((sig, obs_r, p_emp))
    perm_df = pd.DataFrame(rows, columns=['sig_id','obs_r','p_emp']).set_index('sig_id')
    # FDR
    try:
        from statsmodels.stats.multitest import multipletests
        perm_df['fdr'] = multipletests(perm_df['p_emp'], method='fdr_bh')[1]
    except Exception:
        perm_df['fdr'] = np.nan
    out_file = f"{out_prefix}_permutation_top{topN}.tsv"
    perm_df.to_csv(out_file, sep='\t')
    print("[PERM] Wrote:", out_file)
    print("[PERM] elapsed:", time.time() - start)
    return perm_df

def main(args):
    rev = clean_input_rev(args.rev)
    print("[INFO] Loaded reversal file:", args.rev, "rows:", len(rev))
    rev.to_csv(args.out_prefix + ".cleaned.tsv", sep='\t', index=False)
    pert = load_pert_info(args.pert)
    print("[INFO] Loaded pert_info:", args.pert, "rows:", len(pert))

    # Merge on pert_id (left join rev onto pert)
    merged = rev.merge(pert, on='pert_id', how='left', suffixes=('','_pertinfo'))
    # flag unmapped heuristically and if no canonical_smiles/pubchem
    merged['unmapped_flag'] = merged.apply(flag_unmapped, axis=1)
    # also flag if canonical_smiles and pubchem missing
    merged['no_pubchem'] = merged['pubchem_cid'].isna() | (merged['pubchem_cid'].astype(str).str.strip()=='-666') | (merged['pubchem_cid'].astype(str).str.strip()=='')
    merged['no_smiles'] = merged.get('canonical_smiles', pd.Series(['']*len(merged))).astype(str).str.strip()==''

    # final unmapped: either unmapped_flag OR no pubchem & no smiles
    merged['final_unmapped'] = merged['unmapped_flag'] | (merged['no_pubchem'] & merged['no_smiles'])

    # Save merged
    merged.to_csv(args.out_prefix + ".merged.tsv", sep='\t', index=False)

    # Write unmapped prioritized: sort by abs(reversal_score) and by n_signatures
    unmapped = merged[merged['final_unmapped']].copy()
    unmapped['abs_score'] = unmapped['reversal_score'].abs()
    unmapped = unmapped.sort_values(['abs_score','n_signatures'], ascending=[False, False])
    unmapped[['pert_id','pert_iname','reversal_score','n_signatures']].to_csv(
        args.out_prefix + ".unmapped_pertids.tsv", sep='\t', index=False
    )
    print("[INFO] Wrote unmapped list:", args.out_prefix + ".unmapped_pertids.tsv  (count={})".format(len(unmapped)))

    # Aggregate per pert_iname
    agg = aggregate_by_name(merged)
    agg = agg.sort_values(['median_score','total_signatures'], ascending=[False, False])
    agg.to_csv(args.out_prefix + ".aggregated_by_name.tsv", sep='\t', index=False)
    print("[INFO] Wrote aggregated_by_name:", args.out_prefix + ".aggregated_by_name.tsv")

    # Bias adjustments
    agg_adj_log = adjust_for_library_size(agg, method='log')
    agg_adj_power = adjust_for_library_size(agg, method='power')
    # Merge residualized
    agg_resid = residualize_scores(agg)
    # combine into one frame
    comb = agg.copy()
    comb = comb.merge(agg_adj_log[['pert_iname','adj_score']], on='pert_iname', how='left', suffixes=('','_adj_log'))
    comb = comb.merge(agg_resid[['pert_iname','score_resid','onc_flag']], on='pert_iname', how='left')
    comb = comb.merge(agg_adj_power[['pert_iname','adj_score']].rename(columns={'adj_score':'adj_score_power'}), on='pert_iname', how='left')

    # attach representative smiles/pubchem if available from merged (take first non-null)
    rep = merged.groupby('pert_iname').agg({
        'canonical_smiles': lambda x: next((v for v in x if pd.notna(v) and str(v).strip()!='' and str(v)!='-666'), ''),
        'pubchem_cid': lambda x: next((v for v in x if pd.notna(v) and str(v).strip()!='' and str(v)!='-666'), '')
    }).reset_index()
    comb = comb.merge(rep, on='pert_iname', how='left')

    # compute CNS priority & final_score (uses adj_score if present)
    comb = compute_cns_priority(comb)
    comb = comb.sort_values('final_score', ascending=False)
    comb.to_csv(args.out_prefix + ".final_prioritized.tsv", sep='\t', index=False)
    comb.to_csv(args.out_prefix + ".bias_adjusted.tsv", sep='\t', index=False)  # same file for convenience
    print("[INFO] Wrote final prioritized:", args.out_prefix + ".final_prioritized.tsv")

    # quick diagnostics print
    top200 = comb.head(200)
    oncology_frac_top200 = top200['pert_iname'].apply(lambda x: 1 if is_probably_oncology(x) else 0).mean()
    print(f"[DIAG] Top200 oncology-like fraction (heuristic): {oncology_frac_top200:.3f}")

    # Optional permutation testing if provided
    if args.permute and args.permute_drugmat and args.permute_c4:
        print("[INFO] Running permutation testing (this is slow)...")
        do_permutation_test(args.permute_drugmat, args.permute_c4, topN=args.topN_permute, n_perm=args.n_perm,
                            out_prefix=args.out_prefix)
    print("[DONE] All outputs written. Inspect the *.final_prioritized.tsv and unmapped list for manual curation.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--rev', required=True, help='C4A_drug_reversals.tsv path (TSV)')
    p.add_argument('--pert', default='data/raw/GSE92742_Broad_LINCS_pert_info.txt', help='pert_info path (TSV or CSV)')
    p.add_argument('--out_prefix', default='C4A_drug_reversals.optimized', help='output prefix')
    # permutation optional
    p.add_argument('--permute', action='store_true', help='enable permutation testing (requires --permute_drugmat and --permute_c4)')
    p.add_argument('--permute_drugmat', help='drug_gene_matrix.tsv (sig x gene) for permutation')
    p.add_argument('--permute_c4', help='c4a_signature.tsv (gene \t value) for permutation')
    p.add_argument('--topN_permute', type=int, default=200, help='top N signatures for permutation')
    p.add_argument('--n_perm', type=int, default=1000, help='number of permutations for permutation testing')
    args = p.parse_args()
    # Validate permutation args
    if args.permute:
        if not args.permute_drugmat or not args.permute_c4:
            print("[ERROR] --permute requires --permute_drugmat and --permute_c4", file=sys.stderr)
            sys.exit(1)
    main(args)
