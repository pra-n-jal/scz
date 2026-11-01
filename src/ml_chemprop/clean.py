import argparse
import os
import sys
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Clean and isolate top drug-like leads following Chemprop prediction.")
    parser.add_argument("--input", "-i", default="results/predictions_filtered.tsv",
                        help="Input filtered predictions TSV (default: results/predictions_filtered.tsv)")
    parser.add_argument("--out", "-o", default="results/final_shortlist_chemprop.tsv",
                        help="Output shortlist TSV (default: results/final_shortlist_chemprop.tsv)")
    parser.add_argument("--max-mw", type=float, default=600.0,
                        help="Maximum molecular weight (default: 600)")
    parser.add_argument("--min-qed", type=float, default=0.25,
                        help="Minimum quantitative estimate of drug-likeness (default: 0.25)")
    args = parser.parse_args()

    input_path = args.input
    if not os.path.exists(input_path):
        fallbacks = ["Reversal/predictions_filtered.tsv", "predictions_filtered.tsv"]
        for fb in fallbacks:
            if os.path.exists(fb):
                input_path = fb
                break
        else:
            print(f"[ERROR] Input file not found: {args.input}")
            sys.exit(1)

    print(f"[INFO] Reading predictions from {input_path}")
    df = pd.read_csv(input_path, sep="\t")

    mask = pd.Series(True, index=df.index)
    if "is_pains" in df.columns:
        mask &= (df["is_pains"] == False) | (df["is_pains"] == 0) | df["is_pains"].isna()
    if "has_nitro" in df.columns:
        mask &= (df["has_nitro"] == False) | (df["has_nitro"] == 0) | df["has_nitro"].isna()
    if "mw" in df.columns:
        mask &= pd.to_numeric(df["mw"], errors="coerce") < args.max_mw
    if "qed" in df.columns:
        mask &= pd.to_numeric(df["qed"], errors="coerce") > args.min_qed

    keep = df[mask].copy()

    sort_col = "label" if "label" in keep.columns else ("score" if "score" in keep.columns else keep.columns[1])
    if sort_col in keep.columns:
        keep = keep.sort_values(sort_col, ascending=False)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    keep.to_csv(out_path, sep="\t", index=False)
    print(f"[SUCCESS] Filtered {len(df)} compounds -> {len(keep)} drug-like candidates saved to {out_path}")


if __name__ == "__main__":
    main()
