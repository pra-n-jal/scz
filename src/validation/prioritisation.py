import pandas as pd

# load your results
res = pd.read_csv(r"C:\Users\pranj\scz_connectivity_multi\Final_results\Manual\Training\Reversal\C4A_drug_reversals.tsv", sep="\t")

# load mapping (download pert_info from LINCS/CLUE or GEO)
pert = pd.read_csv(r"C:\Users\pranj\scz_connectivity_multi\Final_results\Manual\Training\GSE92742_Broad_LINCS_pert_info.csv", sep="\t", dtype=str)

# merge
m = res.merge(pert, on="pert_id", how="left")

# basic filters and sorting
m['n_signatures'] = pd.to_numeric(m['n_signatures'], errors='coerce').fillna(0).astype(int)
m = m[m['n_signatures'] >= 3]   # robustness threshold
# flag whether mapped to known drugs
m['mapped'] = ~m['pert_iname_y'].isna()
# simple prioritization score (example)
m['priority'] = (m['reversal_score'].astype(float).rank(ascending=False) * 0.6 +
                 m['n_signatures'].rank(ascending=False) * 0.3 +
                 m['mapped'].astype(int) * 0.1)

m = m.sort_values('priority', ascending=False)
m[['pert_id','pert_iname_x','pert_iname_y','moa','reversal_score','single_gene_value','n_signatures']].head(50).to_csv("mapped_prioritized.tsv", sep="\t", index=False)
print("Wrote mapped_prioritized.tsv")
