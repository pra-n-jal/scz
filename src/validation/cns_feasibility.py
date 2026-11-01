from pathlib import Path
import argparse
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def main():
    parser = argparse.ArgumentParser(description="Evaluate Central Nervous System (CNS) MPO and Blood-Brain Barrier feasibility.")
    parser.add_argument("--input", "--mapped", dest="mapped", required=True,
                        help="Input TSV with compound structures or SMILES")
    parser.add_argument("--out", "-o", default="results/final_results_bias_adjusted.tsv",
                        help="Output feasibility TSV")
    args = parser.parse_args()

    input_path = Path(args.mapped)
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {args.mapped}")
        sys.exit(1)

    print(f"[INFO] Reading compound data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    # Detect smiles column
    smiles_col = None
    for c in ["canonical_smiles", "smiles", "SMILES", "can_smiles"]:
        if c in df.columns:
            smiles_col = c
            break

    # Detect ID / name columns
    pert_id_col = "pert_id" if "pert_id" in df.columns else df.columns[0]
    pert_iname_col = "pert_iname" if "pert_iname" in df.columns else ("Compound" if "Compound" in df.columns else pert_id_col)

    rows = []
    for _, r in df.iterrows():
        smi = str(r.get(smiles_col, "")) if smiles_col else ""
        mol = Chem.MolFromSmiles(smi) if smi and smi not in ("nan", "-666", "None", "") else None
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            rotb = Descriptors.NumRotatableBonds(mol)
            lipinski = int((mw < 500) & (logp < 5) & (tpsa < 90))
            cns = int((mw < 450) & (tpsa < 90) & (logp > -0.5) & (logp < 5))
        else:
            mw = logp = tpsa = rotb = lipinski = cns = None

        row_dict = {
            "pert_id": r.get(pert_id_col, ""),
            "pert_iname": r.get(pert_iname_col, ""),
            "MW": mw,
            "LogP": logp,
            "TPSA": tpsa,
            "RotBonds": rotb,
            "Lipinski_OK": lipinski,
            "CNS_pred": cns
        }
        # Carry over original score columns if present
        for sc in ["reversal_score", "spearman", "adj_score", "median_score", "n_signatures"]:
            if sc in r:
                row_dict[sc] = r[sc]
        rows.append(row_dict)

    out_df = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Evaluated CNS parameters for {len(out_df)} compounds -> {out_path}")


if __name__ == "__main__":
    main()
