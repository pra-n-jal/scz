import argparse
import os
import sys
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Summarize directional consistency of drug reversal scores.")
    parser.add_argument("--input", "-i", default="results/FINAL_OUTPUT.tsv",
                        help="Input reversal TSV file (default: results/FINAL_OUTPUT.tsv)")
    parser.add_argument("--out", "-o", default="results/reversal_summary_stats.tsv",
                        help="Output summary TSV file (default: results/reversal_summary_stats.tsv)")
    args = parser.parse_args()

    input_path = args.input
    if not os.path.exists(input_path):
        fallbacks = ["results/C4A_drug_reversals_with_significance.tsv", "C4A_drug_reversals.tsv"]
        for fb in fallbacks:
            if os.path.exists(fb):
                input_path = fb
                break
        else:
            print(f"[ERROR] Could not find input file: {args.input}")
            sys.exit(1)

    print(f"[INFO] Reading {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    score_col = "reversal_score" if "reversal_score" in df.columns else "spearman"
    iname_col = "pert_iname" if "pert_iname" in df.columns else ("Compound" if "Compound" in df.columns else df.columns[0])

    summary = df.groupby(iname_col)[score_col].agg(["mean", "std", "count"]).reset_index()
    summary = summary.rename(columns={"count": "n_observations"})

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved summary statistics for {len(summary)} compounds -> {out_path}")


if __name__ == "__main__":
    main()
