import sys
import argparse
import pandas as pd
from pathlib import Path

DEFAULT_DATASETS = ["GSE21138", "GSE38485", "GSE38484", "GSE48072"]
TOP_N = 5000


def process_datasets(base_dir: Path, datasets=None, top_n=TOP_N):
    if datasets is None:
        datasets = DEFAULT_DATASETS

    for dataset in datasets:
        ds_dir = base_dir / dataset
        de_file = ds_dir / "differential_expression.tsv"

        if not de_file.exists():
            print(f"[DEBUG] {dataset}: {de_file} not found, skipping.")
            continue

        print(f"[INFO] Processing {dataset}...")
        try:
            de_df = pd.read_csv(de_file, sep="\t")
        except Exception as e:
            print(f"[ERROR] Could not read {de_file}: {e}")
            continue

        required_cols = ["gene_symbol", "logFC", "p_value"]
        if not all(col in de_df.columns for col in required_cols):
            print(f"[WARN] {dataset}: differential_expression.tsv missing required columns {required_cols}, skipping.")
            continue

        de_df["ranking_score"] = de_df["logFC"].abs()
        top_genes = de_df.sort_values("ranking_score", ascending=False).head(top_n)

        top_genes_file = ds_dir / "top_genes.tsv"
        top_genes.to_csv(top_genes_file, sep="\t", index=False)
        print(f"[INFO] Saved {top_genes_file}")

        top_genes_mapped = top_genes.groupby("gene_symbol", as_index=False).agg({
            "logFC": "max",
            "p_value": "min",
            "ranking_score": "max"
        })
        top_genes_mapped_file = ds_dir / "top_genes_mapped.tsv"
        top_genes_mapped.to_csv(top_genes_mapped_file, sep="\t", index=False)
        print(f"[INFO] Saved {top_genes_mapped_file}")


def main():
    parser = argparse.ArgumentParser(description="Map and rank top genes from differential expression tables.")
    parser.add_argument("--base-dir", "-b", default="data/raw", help="Base directory containing dataset folders")
    parser.add_argument("--top-n", type=int, default=TOP_N, help="Number of top genes to retain")
    args = parser.parse_args()

    base_p = Path(args.base_dir)
    process_datasets(base_p, top_n=args.top_n)


if __name__ == "__main__":
    main()
