# summarize_signature_consistency.py
import pandas as pd

df = pd.read_csv("C4A_drug_reversals.tsv", sep="\t")
summary = df.groupby("pert_iname")["reversal_score"].agg(["mean", "std", "count"]).reset_index()
summary.to_csv("reversal_summary_stats.tsv", sep="\t", index=False)
print("[INFO] Saved reversal_summary_stats.tsv")
