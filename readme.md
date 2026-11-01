# C4A Schizophrenia Drug Repurposing Pipeline
Finding small molecules that flip the C4A transcriptomic axis to rescue synaptic loss in schizophrenia.



## Running it locally

### Setup
Clone the repo and create a clean environment:

```bash
git clone https://github.com/<your-username>/C4A_Schizophrenia_Repurposing_Project.git
cd C4A_Schizophrenia_Repurposing_Project

conda create -n c4a_pipeline python=3.10 -y
conda activate c4a_pipeline
pip install -r requirements.txt
```

*(Note for Windows users: if pip has trouble installing RDKit wheels, just run `conda install -c conda-forge rdkit`.)*

### Run the tests
You can run the test suite anytime to check that the statistical routines, parsers, and data files are working:
```bash
python -m unittest tests/test_pipeline.py
```

### Running the pipeline step by step

**1. Build the C4A disease signature:**
```bash
python src/signature/build_c_4_a_complement_signature.py \
    --input data/intermediate/C4A_genes_combined.tsv \
    --out data/intermediate/C4A_complement_signature.tsv

# Convert symbols to Entrez IDs
python src/signature/convert_signature_to_entr.py \
    --input data/intermediate/C4A_complement_signature.tsv \
    --out data/intermediate/C4A_signature_entrez.tsv
```

**2. Query LINCS L1000 & calculate drug reversals:**
```bash
python src/lincs/compute_connectivity_lincs.py \
    --sig data/intermediate/C4A_complement_signature.tsv \
    --gctx data/raw/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx \
    --pert_info data/raw/GSE92742_Broad_LINCS_pert_info.txt \
    --out results/connectivity_by_pert.tsv

# Aggregate scores across cell lines to drug level
python src/lincs/C4A_drug_reversals.py \
    --pert-scores results/connectivity_by_pert.tsv \
    --pert-info data/raw/GSE92742_Broad_LINCS_pert_info.txt \
    --out results/C4A_drug_reversals.tsv
```

**3. Chemprop GNN modeling & chemical filtering:**
```bash
# Prepare dataset splits
python src/ml_chemprop/prep_and_split_for_chemprop.py \
    --input data/intermediate/C4A_training_dataset_filled.tsv \
    --out-prefix data/intermediate/chemprop \
    --seed 42

# Filter predictions (strip PAINS, nitro groups, enforce QED > 0.25)
python src/ml_chemprop/filter_and_cluster.py \
    --input results/chemprop_C4A_predictions.csv \
    --out results/predictions_filtered.tsv

# Clean and isolate top candidates
python src/ml_chemprop/clean.py \
    --input results/predictions_filtered.tsv \
    --out results/final_shortlist_chemprop.tsv
```

**4. Permutation FDR testing & CNS scoring:**
```bash
# 10,000-iteration permutation null test
python src/validation/compute_permutation_fdr.py \
    --reversals results/FINAL_OUTPUT.tsv \
    --n_perms 10000 \
    --out results/C4A_with_fdr.tsv

# Score CNS feasibility and BBB permeability
python src/validation/cns_feasibility.py \
    --input results/FINAL_OUTPUT.tsv \
    --out results/final_results_bias_adjusted.tsv
```

---