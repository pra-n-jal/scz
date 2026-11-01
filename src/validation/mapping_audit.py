# mapping_audit.py
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--reversals", required=True)
parser.add_argument("--pert_info", required=True)
parser.add_argument("--out", default="mapping_audit.tsv")
args = parser.parse_args()

rev = pd.read_csv(args.reversals, sep="\t")
info = pd.read_csv(args.pert_info)

merged = rev.merge(info, on="pert_id", how="left")
mapped = merged[merged["canonical_smiles"].notna()]
unmapped = merged[merged["canonical_smiles"].isna()]

print(f"[INFO] Mapped: {len(mapped)} | Unmapped: {len(unmapped)}")
mapped.to_csv("mapped_" + args.out, sep="\t", index=False)
unmapped.to_csv("unmapped_" + args.out, sep="\t", index=False)
merged.to_csv(args.out, sep="\t", index=False)
