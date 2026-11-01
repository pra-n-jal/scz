from pathlib import Path
import pandas as pd
import numpy as np
import logging
import argparse

DATASETS = ["GSE38484", "GSE38485", "GSE21138", "GSE48072"]
MANUAL_DIR = Path(".")
DIFF_FILENAME = "differential_expression.tsv"
OUTPUT_SIGNATURE = Path("C4A_complement_signature.tsv")
EXTRACT_DIR = Path("extracted")

PREDEFINED_COMPLEMENT_GENES = [
    "C1QA", "C1QB", "C1QC",
    "C2", "C3", "C4A", "C4B",
    "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
    "CFH", "CFB", "CR1",
]

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

POSSIBLE_GENE_COLS = ["gene", "Gene", "gene_symbol", "GeneSymbol", "GENE", "symbol", "hgnc_symbol"]
POSSIBLE_LOGFC_COLS = ["logFC", "log2FoldChange", "log2Fold_Change", "log2fc", "logFC_baseline", "mean_logFC", "combined_z", "score", "value"]
POSSIBLE_P_COLS = ["P.Value", "pval", "p_value", "p.value", "adj.P.Val", "padj", "FDR", "qval"]


def detect_columns(df):
    gene_col = None
    for c in POSSIBLE_GENE_COLS:
        if c in df.columns:
            gene_col = c
            break
    if gene_col is None:
        for c in df.columns:
            sample = df[c].astype(str).head(50).tolist()
            alpha_frac = sum(any(ch.isalpha() for ch in s) for s in sample) / len(sample) if sample else 0
            if alpha_frac > 0.6:
                gene_col = c
                break

    logfc_col = None
    for c in POSSIBLE_LOGFC_COLS:
        if c in df.columns:
            logfc_col = c
            break
    if logfc_col is None:
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        for c in numeric_cols:
            if c.lower() not in [x.lower() for x in POSSIBLE_P_COLS]:
                logfc_col = c
                break

    p_col = None
    for c in POSSIBLE_P_COLS:
        if c in df.columns:
            p_col = c
            break
    return gene_col, logfc_col, p_col


def read_diff_file(path: Path):
    if not path.exists():
        logging.warning(f"File not found: {path}")
        return None
    try:
        df = pd.read_csv(path, sep="\t", index_col=False)
    except Exception as e:
        logging.error(f"Could not read {path}: {e}")
        return None
    gene_col, logfc_col, p_col = detect_columns(df)
    if gene_col is None or logfc_col is None:
        logging.error(f"Could not detect gene or logFC columns in {path}. Columns: {list(df.columns)}")
        return None
    df2 = df[[gene_col, logfc_col]].copy()
    df2.columns = ["gene", "logFC"]
    df2["gene"] = df2["gene"].astype(str).str.upper().str.strip()
    df2["logFC"] = pd.to_numeric(df2["logFC"], errors="coerce")
    return df2.dropna(subset=["gene"])


def build_signature_from_input_file(input_path: Path, complement_genes=None):
    if not input_path.exists():
        logging.error(f"Input file not found: {input_path}")
        return None
    try:
        df = pd.read_csv(input_path, sep="\t", index_col=False)
    except Exception:
        df = pd.read_csv(input_path, index_col=False)

    gene_col, val_col, _ = detect_columns(df)
    if gene_col is None or val_col is None:
        logging.error(f"Could not detect gene or value column in {input_path}")
        return None

    df_out = df.copy()
    df_out["gene"] = df_out[gene_col].astype(str).str.upper().str.strip()
    if complement_genes:
        complement_genes_up = [g.upper() for g in complement_genes]
        matching = df_out[df_out["gene"].isin(complement_genes_up)]
        if len(matching) > 0:
            df_out = matching

    # If already has combined_z, keep it
    if "combined_z" in df_out.columns:
        df_out["combined_z"] = pd.to_numeric(df_out["combined_z"], errors="coerce")
    elif val_col:
        df_out["combined_z"] = pd.to_numeric(df_out[val_col], errors="coerce")
    else:
        df_out["combined_z"] = 1.0

    if "mean_logFC" not in df_out.columns:
        df_out["mean_logFC"] = pd.to_numeric(df_out[val_col], errors="coerce")

    df_out = df_out.sort_values("combined_z", ascending=False).drop_duplicates(subset=["gene"])
    cols = ["gene", "mean_logFC", "combined_z"] + [c for c in df_out.columns if c not in ("gene", "mean_logFC", "combined_z")]
    return df_out[cols].set_index("gene")


def build_signature(datasets, complement_genes):
    complement_genes_up = [g.upper() for g in complement_genes]
    EXTRACT_DIR.mkdir(exist_ok=True)
    per_ds = {}
    for ds in datasets:
        ds_path = MANUAL_DIR / ds
        diff_path = ds_path / DIFF_FILENAME
        df = read_diff_file(diff_path)
        if df is None:
            logging.warning(f"Skipping dataset {ds} (no usable diff file)")
            continue

        df_sub = df[df['gene'].isin(complement_genes_up)].set_index('gene')
        outp = EXTRACT_DIR / f"{ds}_complement_extracted.tsv"
        df_sub.to_csv(outp, sep="\t")
        logging.info(f"Extracted {len(df_sub)} complement genes from {ds} -> {outp}")
        per_ds[ds] = df_sub['logFC']

    if not per_ds:
        logging.error("No datasets yielded complement genes.")
        return None

    all_genes = sorted({g for s in per_ds.values() for g in s.index})
    mat = pd.DataFrame(index=all_genes, columns=per_ds.keys(), dtype=float)
    for ds, series in per_ds.items():
        mat.loc[series.index, ds] = series

    zmat = mat.copy()
    for ds in zmat.columns:
        col = zmat[ds]
        mean = col.mean(skipna=True)
        std = col.std(skipna=True)
        if pd.isna(std) or std == 0:
            zmat[ds] = (col - mean)
        else:
            zmat[ds] = (col - mean) / std

    z_sum = zmat.sum(axis=1, skipna=True)
    k_nonmiss = zmat.notna().sum(axis=1).astype(float)
    combined_z = z_sum / np.sqrt(k_nonmiss)

    mean_logfc = mat.mean(axis=1, skipna=True)

    result = pd.DataFrame({
        'gene': mat.index,
        'mean_logFC': mean_logfc.values,
        'combined_z': combined_z.values,
        'n_datasets': k_nonmiss.values
    }).set_index('gene')

    for ds in mat.columns:
        result[f'logFC_{ds}'] = mat[ds]

    result = result.sort_values('combined_z', ascending=False)
    return result


def main(args):
    if args.genes_file:
        if not Path(args.genes_file).exists():
            logging.error(f"Genes file not found: {args.genes_file}")
            return
        user_genes = [x.strip() for x in Path(args.genes_file).read_text().splitlines() if x.strip()]
        complement_genes = user_genes
        logging.info(f"Loaded {len(complement_genes)} complement genes from {args.genes_file}")
    else:
        complement_genes = PREDEFINED_COMPLEMENT_GENES
        logging.info(f"Using built-in complement gene list ({len(complement_genes)} genes)")

    output_path = Path(args.output if args.output else args.out)

    if args.input:
        input_path = Path(args.input)
        logging.info(f"Building signature from input file: {input_path}")
        signature = build_signature_from_input_file(input_path, complement_genes)
    else:
        datasets = args.datasets if args.datasets else DATASETS
        logging.info(f"Datasets: {datasets}")
        signature = build_signature(datasets, complement_genes)

    if signature is None:
        logging.error("Signature build failed.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    signature.to_csv(output_path, sep="\t")
    logging.info(f"Wrote signature to {output_path} (rows={len(signature)})")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Build combined C4A / complement signature from differential_expression.tsv files or combined table")
    ap.add_argument('--input', '-i', help='Optional path to combined differential expression table', default=None)
    ap.add_argument('--datasets', nargs='+', help='List of dataset folder names to include (overrides config)', default=None)
    ap.add_argument('--genes-file', help='Optional file (one gene per line) with complement genes to use', default=None)
    ap.add_argument('--output', '--out', dest='output', help='Output signature TSV', default=str(OUTPUT_SIGNATURE))
    args = ap.parse_args()
    args.out = args.output
    main(args)
