import pandas as pd
df = pd.read_csv("Reversal/predictions_filtered.tsv", sep='—t')
keep = df[(df['is_pains']==False) & (df['has_nitro']==False) & (df['mw']<600) & (df['qed']>0.25)]
keep.sort_values('label', ascending=False).to_csv("Reversal/predictions_filtered_top.tsv", sep='—t', index=False)
print("Kept:", len(keep))
