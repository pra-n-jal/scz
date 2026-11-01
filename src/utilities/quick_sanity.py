import argparse
import sys
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Quick sanity check mapping reversal results to LINCS perturbagen metadata.")
    parser.add_argument('--reversal', default='results/C4A_drug_reversals_with_significance.tsv',
                        help='Reversal TSV file')
    parser.add_argument('--pert-info', default='data/raw/GSE92742_Broad_LINCS_pert_info.txt',
                        help='LINCS pert_info file (TSV)')
    parser.add_argument('--out', default='results/top20_annotated.tsv',
                        help='Output annotated top 20 TSV')
    args = parser.parse_args()

    rev_path = Path(args.reversal)
    if not rev_path.exists():
        alt = Path("results/C4A_drug_reversals.tsv")
        if alt.exists():
            rev_path = alt
        else:
            print(f"[ERROR] Reversal file not found: {args.reversal}")
            sys.exit(1)

    pert_path = Path(args.pert_info)
    if not pert_path.exists():
        print(f"[ERROR] pert_info file not found: {args.pert_info}")
        sys.exit(1)

    rev = pd.read_csv(rev_path, sep='\t')
    pert_info = pd.read_csv(pert_path, sep='\t', low_memory=False)

    print('Total rows in reversal:', len(rev))
    iname_col = 'pert_iname' if 'pert_iname' in rev.columns else rev.columns[0]
    print(f'Unique {iname_col}:', rev[iname_col].nunique())

    merged = rev.merge(
        pert_info[['pert_iname', 'pert_type', 'moa', 'pert_id']],
        left_on=iname_col,
        right_on='pert_iname',
        how='left'
    )

    mapped = merged['pert_type'].notnull().sum()
    print(f'Mapped to pert_info: {mapped} / {len(merged)} ({mapped / len(merged):.1%})')
    print('pert_type counts:\n', merged['pert_type'].value_counts(dropna=False).head(20))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sort_col = 'reversal_score' if 'reversal_score' in merged.columns else ('spearman' if 'spearman' in merged.columns else merged.columns[1])
    merged.sort_values(sort_col).head(20).to_csv(out, sep='\t', index=False)
    print('Wrote', out)


if __name__ == '__main__':
    main()
