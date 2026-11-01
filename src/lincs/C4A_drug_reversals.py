from pathlib import Path
import argparse
import sys
import os
import pandas as pd


def to_numeric_if_exists(df, col):
    if col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser(description="Aggregate LINCS connectivity scores into drug-level reversal tables.")
    parser.add_argument('--pert-scores', default='results/connectivity_by_pert.tsv',
                        help='Input per-signature or per-pert scores TSV')
    parser.add_argument('--col-meta', default='lincs_col_meta_parsed.tsv',
                        help='Parsed LINCS column metadata TSV (optional if pert_id is already present)')
    parser.add_argument('--pert-info', default='data/raw/GSE92742_Broad_LINCS_pert_info.txt',
                        help='Broad LINCS perturbagen metadata (TSV)')
    parser.add_argument('--out', default='results/C4A_drug_reversals.tsv',
                        help='Output drug-level reversal TSV')
    args = parser.parse_args()

    if not os.path.exists(args.pert_scores):
        alt = "Reversal/connectivity_by_pert.tsv"
        if os.path.exists(alt):
            args.pert_scores = alt
        else:
            print(f"[FATAL] Perturbation scores file not found: {args.pert_scores}")
            sys.exit(1)

    print(f'[INFO] Loading perturbation scores from {args.pert_scores}...')
    pert_scores = pd.read_csv(args.pert_scores, sep='\t', dtype=str, low_memory=False)
    to_numeric_if_exists(pert_scores, 'spearman')
    to_numeric_if_exists(pert_scores, 'reversal_score')
    pert_scores.columns = [c.strip() for c in pert_scores.columns]

    col_meta = None
    if os.path.exists(args.col_meta):
        col_meta = pd.read_csv(args.col_meta, sep='\t', dtype=str, low_memory=False)
        col_meta.columns = [c.strip() for c in col_meta.columns]
    else:
        print(f"[DEBUG] Optional col_meta '{args.col_meta}' not found; skipping sig_id fallback parsing.")

    pert_info = None
    if os.path.exists(args.pert_info):
        pert_info = pd.read_csv(args.pert_info, sep='\t', dtype=str, low_memory=False)
        pert_info.columns = [c.strip() for c in pert_info.columns]
    else:
        print(f"[DEBUG] Optional pert_info '{args.pert_info}' not found.")

    if 'pert_id' not in pert_scores.columns:
        if 'sig_id' in pert_scores.columns and col_meta is not None:
            print("[DEBUG] 'pert_id' not found in pert_scores — attempting to map via col_meta using sig_id...")
            merged_tmp = pert_scores.merge(col_meta, left_on='sig_id', right_on='id', how='left')
            if 'full_pert_batch' in merged_tmp.columns and merged_tmp['full_pert_batch'].notna().any():
                merged_tmp['pert_id'] = merged_tmp['full_pert_batch']
                pert_scores = merged_tmp
                print('[DEBUG] Derived pert_id from full_pert_batch in col_meta.')
            elif 'pert_id' in merged_tmp.columns and merged_tmp['pert_id'].notna().any():
                pert_scores = merged_tmp
                print("[DEBUG] Derived pert_id from col_meta 'pert_id' column.")
            else:
                print('[FATAL] Could not derive pert_id from col_meta.')
                sys.exit(1)
        else:
            if 'pert_iname' in pert_scores.columns:
                pert_scores['pert_id'] = pert_scores['pert_iname']
            else:
                print("[FATAL] pert_scores lacks 'pert_id', 'sig_id', and 'pert_iname'. Cannot proceed.")
                sys.exit(1)

    pert_scores['pert_id'] = pert_scores['pert_id'].astype(str).fillna('')

    # If pert_info is available, map compound annotations
    if pert_info is not None and 'pert_id' in pert_info.columns:
        pert_info_cp = pert_info[pert_info.get('pert_type', '') == 'trt_cp'].copy() if 'pert_type' in pert_info.columns else pert_info.copy()
        if pert_info_cp.empty:
            pert_info_cp = pert_info.copy()
        print(f'[DEBUG] Merging pert_scores with pert_info on pert_id...')
        merged = pert_scores.merge(pert_info_cp, on='pert_id', how='left')
    else:
        merged = pert_scores.copy()

    # Determine pert_iname
    iname_col = None
    for cand in ['pert_iname_y', 'pert_iname_x', 'pert_iname', 'compound', 'Compound']:
        if cand in merged.columns and merged[cand].notna().any():
            iname_col = cand
            break

    if iname_col:
        merged['pert_iname'] = merged[iname_col]
    elif 'pert_id' in merged.columns:
        merged['pert_iname'] = merged['pert_id']

    compound_data = merged.dropna(subset=['pert_iname']).copy()
    if compound_data.empty:
        print('[ERROR] No perturbation rows with compound names present. Exiting.')
        sys.exit(1)

    agg_cols = {}
    if 'spearman' in compound_data.columns:
        to_numeric_if_exists(compound_data, 'spearman')
        agg_cols['spearman'] = 'mean'
    if 'reversal_score' in compound_data.columns:
        to_numeric_if_exists(compound_data, 'reversal_score')
        agg_cols['reversal_score'] = 'mean'
    if 'single_gene_value' in compound_data.columns and compound_data['single_gene_value'].notna().any():
        agg_cols['single_gene_value'] = 'first'

    group_keys = ['pert_id', 'pert_iname'] if 'pert_id' in compound_data.columns else ['pert_iname']
    agg = compound_data.groupby(group_keys).agg(agg_cols if agg_cols else {group_keys[0]: 'size'}).reset_index()
    agg['n_signatures'] = compound_data.groupby(group_keys).size().values

    if 'reversal_score' in agg.columns:
        agg['reversal_score'] = pd.to_numeric(agg['reversal_score'], errors='coerce')
        agg = agg.sort_values('reversal_score', ascending=False)
    elif 'spearman' in agg.columns:
        agg['spearman'] = pd.to_numeric(agg['spearman'], errors='coerce')
        agg = agg.sort_values('spearman', ascending=False)

    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(outpath, sep='\t', index=False)
    print(f'[SUCCESS] Wrote pert-level ranked reversal table -> {outpath}')
    print(f'[INFO] Total unique perturbagens ranked: {len(agg)}')


if __name__ == '__main__':
    main()
