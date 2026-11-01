import argparse
import os
import sys
from pathlib import Path
import pandas as pd
from mygene import MyGeneInfo


def parse_args():
    parser = argparse.ArgumentParser(description="Map gene symbols in signature to Entrez IDs via MyGene.")
    parser.add_argument("--input", "-i", default="data/intermediate/C4A_complement_signature.tsv",
                        help="Input signature TSV with gene symbols (default: data/intermediate/C4A_complement_signature.tsv)")
    parser.add_argument("--out", "-o", default="data/intermediate/C4A_signature_entrez.tsv",
                        help="Output signature TSV with Entrez IDs (default: data/intermediate/C4A_signature_entrez.tsv)")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = args.input

    if not os.path.exists(input_path):
        # Fallback to current directory if not found in data/intermediate
        alt = os.path.basename(input_path)
        if os.path.exists(alt):
            input_path = alt
        else:
            print(f"[ERROR] Signature file not found at '{args.input}'.")
            sys.exit(1)

    print(f"[INFO] Loading signature from {input_path}")
    sig = pd.read_csv(input_path, sep="\t")

    # Detect gene column
    gene_col = None
    for c in sig.columns:
        if c.strip().lower() in ("gene", "gene_symbol", "symbol", "genes"):
            gene_col = c
            break
    if gene_col is None:
        gene_col = sig.columns[0]

    # Detect score/z column
    val_col = None
    for c in sig.columns:
        if c.strip().lower() in ("combined_z", "z_score", "score", "logfc", "mean_logfc", "value"):
            val_col = c
            break
    if val_col is None:
        val_col = sig.columns[1] if len(sig.columns) > 1 else gene_col

    genes = sig[gene_col].astype(str).str.strip().tolist()
    print(f"[INFO] Querying MyGene for {len(genes)} gene symbols...")
    mg = MyGeneInfo()
    try:
        out = mg.querymany(genes, scopes="symbol", fields="entrezgene", species="human")
        mapping = {d["query"]: d.get("entrezgene", None) for d in out}
    except Exception as e:
        print(f"[ERROR] MyGene query failed: {e}")
        sys.exit(1)

    sig["entrez_id"] = sig[gene_col].map(mapping)
    sig = sig.dropna(subset=["entrez_id"])

    out_df = pd.DataFrame({
        "gene": sig["entrez_id"].astype(str),
        "combined_z": sig[val_col]
    })

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved mapped signature with {len(out_df)} genes -> {out_path}")


if __name__ == "__main__":
    main()
