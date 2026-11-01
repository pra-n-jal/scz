import argparse
import os
import sys
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Multi-attribute prioritization based on reversal score and signature support.")
    parser.add_argument("--reversals", default="results/C4A_drug_reversals_with_significance.tsv",
                        help="Reversal score table (TSV)")
    parser.add_argument("--pert-info", default="data/raw/GSE92742_Broad_LINCS_pert_info.txt",
                        help="LINCS perturbagen info table (TSV/CSV)")
    parser.add_argument("--out", default="results/mapped_prioritized.tsv",
                        help="Output path for prioritized candidates")
    parser.add_argument("--min-signatures", type=int, default=3,
                        help="Minimum number of signatures required (default: 3)")
    parser.add_argument("--top-k", type=int, default=50,
                        help="Top K candidates to output (default: 50)")
    args = parser.parse_args()

    if not os.path.exists(args.reversals):
        # Fallback check
        alt = "results/C4A_drug_reversals.tsv"
        if os.path.exists(alt):
            args.reversals = alt
        else:
            print(f"[ERROR] Reversals file not found: {args.reversals}")
            sys.exit(1)

    print(f"[INFO] Loading reversals from {args.reversals}")
    res = pd.read_csv(args.reversals, sep="\t")

    if os.path.exists(args.pert_info):
        print(f"[INFO] Loading perturbagen metadata from {args.pert_info}")
        sep = "\t" if args.pert_info.endswith((".tsv", ".txt")) else ","
        pert = pd.read_csv(args.pert_info, sep=sep, dtype=str, low_memory=False)
        m = res.merge(pert, on="pert_id", how="left")
    else:
        print(f"[WARN] pert_info not found at {args.pert_info}. Proceeding with existing columns.")
        m = res.copy()

    # Determine signature count and reversal score columns
    sig_col = "n_signatures" if "n_signatures" in m.columns else "Signature_Count"
    if sig_col in m.columns:
        m[sig_col] = pd.to_numeric(m[sig_col], errors="coerce").fillna(0).astype(int)
        m = m[m[sig_col] >= args.min_signatures]

    iname_col = "pert_iname_y" if "pert_iname_y" in m.columns else ("pert_iname" if "pert_iname" in m.columns else "Compound")
    m["mapped"] = ~m[iname_col].isna() if iname_col in m.columns else 1

    rev_col = "reversal_score" if "reversal_score" in m.columns else "LINCS_Reversal_Score"
    if rev_col in m.columns:
        m["priority"] = (
            m[rev_col].astype(float).rank(ascending=False) * 0.6 +
            (m[sig_col].rank(ascending=False) * 0.3 if sig_col in m.columns else 0.0) +
            m["mapped"].astype(int) * 0.1
        )
        m = m.sort_values("priority", ascending=False)

    out_p = Path(args.out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    m.head(args.top_k).to_csv(out_p, sep="\t", index=False)
    print(f"[INFO] Wrote {min(args.top_k, len(m))} prioritized compounds -> {out_p}")


if __name__ == "__main__":
    main()
