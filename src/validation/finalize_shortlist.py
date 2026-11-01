#!/usr/bin/env python3
"""
finalize_shortlist.py
- Input: Reversal/predictions_filtered.tsv
- Pert info: GSE92742_Broad_LINCS_pert_info.txt (tab)
- Output: Reversal/final_shortlist.tsv (selected, annotated with pert_id if matched)

Usage:
  python Reversal/finalize_shortlist.py --pred Reversal/predictions_filtered_top.tsv --pert GSE92742_Broad_LINCS_pert_info.txt --out Reversal/final_shortlist.tsv --per_cluster 1 --qed 0.35 --mw 550
"""
import argparse, pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

def inchikey_from_smiles(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None: return None
        return AllChem.MolToInchiKey(m)
    except:
        return None

def load_pred(path):
    df = pd.read_csv(path, sep='\t')
    # ensure expected columns
    for c in ['smiles','label','qed','mw','is_pains','has_nitro','cluster_id']:
        if c not in df.columns:
            raise ValueError(f"Missing column {c} in predictions file")
    return df

def map_to_pert(df_pred, pert_info_path):
    pert = pd.read_csv(pert_info_path, sep='\t', low_memory=False)
    # find a column that looks like smiles or canonical_smiles
    smiles_col = None
    for c in pert.columns:
        if 'smiles' in c.lower() or 'canonical' in c.lower():
            smiles_col = c
            break
    if smiles_col:
        pert['inchikey'] = pert[smiles_col].apply(lambda s: inchikey_from_smiles(str(s)) if pd.notna(s) and str(s).strip() else None)
    else:
        pert['inchikey'] = None

    # make mapping dict from inchikey -> list of pert rows
    pert_map = pert[['pert_id','pert_iname','inchikey']].drop_duplicates()
    return pert_map

def main(args):
    pred = load_pred(args.pred)
    # filters
    pred = pred[(pred['qed'] >= args.qed) & (pred['mw'] <= args.mw)]
    pred = pred[(pred['is_pains'] == False) & (pred['has_nitro'] == False)]
    # compute inchikeys if not present
    if 'inchikey' not in pred.columns:
        pred['inchikey'] = pred['smiles'].apply(inchikey_from_smiles)
    # select top per cluster
    per_cluster = args.per_cluster
    selected = pred.sort_values('label', ascending=False).groupby('cluster_id').head(per_cluster).reset_index(drop=True)
    # map to pert_info
    pert_map = map_to_pert(selected, args.pert)
    merged = selected.merge(pert_map, how='left', on='inchikey')
    # if no match by inchikey, try matching by name (pert_iname) later manually
    out = Path(args.out)
    merged.to_csv(out, sep='\t', index=False)
    print("[INFO] Wrote shortlist to", out)
    print("[INFO] Summary: n_selected=", len(merged))
    print(merged[['smiles','label','qed','mw','cluster_id','pert_id','pert_iname','inchikey']].head(50))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--pred', required=True)
    p.add_argument('--pert', required=True)
    p.add_argument('--out', default='Reversal/final_shortlist.tsv')
    p.add_argument('--per_cluster', type=int, default=1)
    p.add_argument('--qed', type=float, default=0.35)
    p.add_argument('--mw', type=float, default=550)
    args = p.parse_args()
    main(args)
