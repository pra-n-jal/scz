# cns_feasibility.py
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mapped", required=True)
parser.add_argument("--out", default="cns_feasibility.tsv")
args = parser.parse_args()

df = pd.read_csv(args.mapped, sep="\t")

rows = []
for _, r in df.iterrows():
    smi = str(r.get("canonical_smiles", ""))
    mol = Chem.MolFromSmiles(smi)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        lipinski = int((mw < 500) & (logp < 5) & (tpsa < 90))
        cns = int((mw < 450) & (tpsa < 90) & (logp > -0.5) & (logp < 5))
    else:
        mw = logp = tpsa = rotb = lipinski = cns = None
    rows.append({
        "pert_id": r["pert_id"],
        "pert_iname": r["pert_iname"],
        "MW": mw, "LogP": logp, "TPSA": tpsa,
        "RotBonds": rotb, "Lipinski_OK": lipinski, "CNS_pred": cns
    })
out = pd.DataFrame(rows)
out.to_csv(args.out, sep="\t", index=False)
print(f"[INFO] Saved {args.out}")
