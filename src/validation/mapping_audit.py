import argparse
import sys
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Audit compound mapping coverage against LINCS perturbagen metadata.")
    parser.add_argument("--reversals", default="results/FINAL_OUTPUT.tsv",
                        help="Reversal results table (TSV)")
    parser.add_argument("--pert_info", default="data/raw/GSE92742_Broad_LINCS_pert_info.txt",
                        help="LINCS perturbagen info table (TSV or CSV)")
    parser.add_argument("--out", default="results/mapping_audit.tsv",
                        help="Output audit table TSV")
    args = parser.parse_args()

    rev_path = Path(args.reversals)
    if not rev_path.exists():
        print(f"[ERROR] Reversals file not found: {args.reversals}")
        sys.exit(1)

    pert_path = Path(args.pert_info)
    if not pert_path.exists():
        print(f"[ERROR] Pert_info file not found: {args.pert_info}")
        sys.exit(1)

    print(f"[INFO] Auditing {rev_path} against {pert_path}...")
    rev = pd.read_csv(rev_path, sep="\t")
    sep = "\t" if pert_path.suffix in [".txt", ".tsv"] else ","
    info = pd.read_csv(pert_path, sep=sep, low_memory=False)

    merged = rev.merge(info, on="pert_id", how="left")
    smiles_col = "canonical_smiles" if "canonical_smiles" in merged.columns else ("smiles" if "smiles" in merged.columns else None)

    if smiles_col:
        mapped = merged[merged[smiles_col].notna() & (merged[smiles_col] != "-666")]
        unmapped = merged[merged[smiles_col].isna() | (merged[smiles_col] == "-666")]
    else:
        mapped = merged
        unmapped = pd.DataFrame(columns=merged.columns)

    print(f"[INFO] Mapped: {len(mapped)} | Unmapped: {len(unmapped)}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapped_out = out_path.parent / f"mapped_{out_path.name}"
    unmapped_out = out_path.parent / f"unmapped_{out_path.name}"

    mapped.to_csv(mapped_out, sep="\t", index=False)
    unmapped.to_csv(unmapped_out, sep="\t", index=False)
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Saved audit results -> {out_path}")


if __name__ == "__main__":
    main()
