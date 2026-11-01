# Data Directory Structure

This folder contains local data used by the C4A Schizophrenia drug repurposing pipeline. Massive public data files (such as LINCS L1000 GCTX matrices) are excluded from version control via `.gitignore`.

## Directory Overview

* **`data/raw/`**: Raw benchmark tables and downloaded perturbational databases.
  * `GSE92742_Broad_LINCS_pert_info.txt`: LINCS L1000 small molecule perturbagen metadata (available from [NCBI GEO GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)).
  * `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`: Optional Level 5 moderate z-score expression matrix (~1.2 GB compressed).
  * GEO Patient Microarray Folders (`GSE38484/`, `GSE38485/`, `GSE21138/`, `GSE48072/`): Each containing `differential_expression.tsv` tables.
* **`data/intermediate/`**: Intermediate processed tables:
  * `C4A_genes_combined.tsv`: Consolidated patient-cohort differential expression table.
  * `C4A_complement_signature.tsv`: Multi-cohort combined z-score consensus signature.
  * `C4A_signature_entrez.tsv`: Entrez-mapped signature.
* **`data/outputs/`**: Pipeline runs and export artifacts.

## Pre-computed Results
Note that all final results and prioritized compound tables are already generated and checked into the [`results/`](../results) directory.
