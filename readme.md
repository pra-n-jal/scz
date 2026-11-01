# C4A Schizophrenia Drug Repurposing Pipeline
### AI-Driven Inversion of the Complement C4A Transcriptomic Axis & Antipsychotics Discovery

This repository contains the complete, modular research pipeline for computational drug repurposing targeting **Complement Component 4A (C4A)** in schizophrenia. By integrating multi-cohort patient transcriptomics with the Broad Institute LINCS L1000 dataset and message-passing graph neural networks (Chemprop), the pipeline identifies small molecules capable of reversing pathological C4A-driven synaptic pruning and neuroinflammation.

Crucially, **the unbiased computational method organically re-discovered clinically established antipsychotics** (such as *thioridazine* and *nor-1-chlorpromazine*), their downstream intracellular effector kinase (**GSK3β**), neuroleptic extrapyramidal rescue therapeutics (*biperiden*), and psychiatric anxiolytics (*lorazepam*, *alprazolam*), validating the biological integrity of the transcriptomic inversion approach.

---

## Repository Structure

```text
C4A_Schizophrenia_Repurposing_Project/
├── docs/                                          # Scientific papers, methodology & background
│   ├── METHODOLOGY_AND_ANTIPSYCHOTICS_DISCOVERY.md# In-depth report on how the method uncovered antipsychotics
│   ├── Research_Paper_Original.pdf                # Original publication paper
│   ├── Research_In_Depth_Analysis.pdf             # Comprehensive 9-page in-depth research analysis
│   ├── Synaptic_Pruning_SCZ_C4A_Background.pdf    # Biological background on C4A and synaptic elimination
│   └── Visuals.html                               # Interactive browser-based methodology & candidate visualization
├── results/                                       # Compiled findings and numerical output datasets
│   ├── antipsychotics_and_cns_candidates.tsv      # Extracted antipsychotic, psychiatric & CNS benchmark hits
│   ├── FINAL_OUTPUT.tsv                           # Full 893-compound LINCS mapped reversal ranking
│   ├── final_results_bias_adjusted.tsv            # Bias-adjusted, signature-weighted prioritization with BBB/CNS
│   ├── final_shortlist_chemprop.tsv               # Chemprop neural network candidates after PAINS/QED filtering
│   ├── C4A_with_fdr.tsv                           # 10,000-iteration permutation null test p-values and FDR
│   ├── C4A_drug_reversals_with_significance.tsv   # Drug-level reversal table with parametric significance
│   └── chemprop_C4A_predictions.csv               # Raw GNN model predictions across chemical library
├── src/                                           # Compiled, modular Python source codebase
│   ├── signature/                                 # Stage 1: Consensus C4A signature construction
│   │   ├── build_c_4_a_complement_signature.py    # Derives directionally consistent signature across GEO cohorts
│   │   └── convert_signature_to_entr.py           # Translates gene symbols to Entrez IDs for connectivity tools
│   ├── lincs/                                     # Stage 2: LINCS L1000 connectivity & reversal analysis
│   │   ├── extract_lincs_metadata.py              # Parses GCTX column metadata and perturbation annotations
│   │   ├── compute_connectivity_lincs.py          # Calculates non-parametric reversal scores against LINCS
│   │   ├── C4A_drug_reversals.py                  # Aggregates scores across cell lines and signature replicates
│   │   ├── analyse_reversal.py                    # Mechanism-of-action and compound ranking
│   │   ├── gctx.py                                # Low-level HDF5/GCTX matrix extraction utility
│   │   └── map_to_lincs.py                        # InChIKey mapping to LINCS perturbational metadata
│   ├── ml_chemprop/                               # Stage 3: Chemprop Graph Neural Network modeling
│   │   ├── prep_and_split_for_chemprop.py         # Formats molecular graphs and generates training splits
│   │   ├── filter_and_cluster.py                  # Calculates molecular descriptors, PAINS flags & Butina clusters
│   │   ├── clean.py                               # Enforces QED, molecular weight, and nitro group filters
│   │   └── reversal.py                            # Merges GNN predictions with disease tables for final ranking
│   ├── validation/                                # Stage 4: Statistical validation & CNS feasibility
│   │   ├── compute_permutation_fdr.py             # 10,000-iteration vector permutation null & empirical FDR
│   │   ├── compute_reversal_significance.py       # Empirical p-values across compound perturbation sets
│   │   ├── summarize_signature_consistency.py     # Cross-cohort directional consistency metrics
│   │   ├── cns_feasibility.py                     # Central Nervous System MPO and BBB permeability scoring
│   │   ├── prioritisation.py                      # Multi-attribute ranking (reversal + CNS feasibility)
│   │   ├── finalize_shortlist.py                  # Shortlist assembly combining ML, CNS, and filters
│   │   ├── mapping_audit.py                       # Audit report for mapped vs unmapped perturbagens
│   │   └── optimize_c4a_reversals.py              # Optimizes scoring weights against biological benchmarks
│   └── utilities/                                 # Supporting utilities
│       ├── map_genes.py                           # BiomaRt / MyGene gene ID conversion
│       ├── check_genes.py                         # Validates gene presence across microarrays
│       ├── fetch_smiles_pubchem.py                # Resolves chemical names/CIDs to canonical SMILES
│       ├── map_lincs_to_smiles.py                 # Annotates LINCS perturbations with SMILES strings
│       ├── quick_sanity.py                        # Verification check joining reversal tables to LINCS
│       └── manual_drug_lookup.py                  # Direct queries for specific perturbagens
└── requirements.txt                               # Python package dependencies
```

---

## Quick Start & Installation

### 1. Environment Setup

```bash
git clone <repo-url>
cd C4A_Schizophrenia_Repurposing_Project

# Create and activate a clean environment
conda create -n c4a_pipeline python=3.10 -y
conda activate c4a_pipeline

# Install core dependencies
pip install -r requirements.txt
```

*Note on RDKit:* For cheminformatics filtering, if `pip install rdkit` encounters binary issues on Windows, install via Conda:
```bash
conda install -c conda-forge rdkit
```

---

## Scientific Workflow & Execution Guide

### Step 1: Build the Consensus C4A Disease Signature
Derives a robust core signature from five independent GEO patient datasets (post-mortem brain BA46, cortical regions, and whole blood):
```bash
python src/signature/build_c_4_a_complement_signature.py \
    --input data/intermediate/C4A_genes_combined.tsv \
    --out data/intermediate/C4A_complement_signature.tsv
```
To convert symbols to Entrez IDs:
```bash
python src/signature/convert_signature_to_entr.py \
    --input data/intermediate/C4A_complement_signature.tsv \
    --out data/intermediate/C4A_signature_entrez.tsv
```

### Step 2: Query the LINCS L1000 Connectivity Matrix
Calculates non-parametric rank connectivity scores against ~20,000 small-molecule perturbation experiments:
```bash
python src/lincs/compute_connectivity_lincs.py \
    --sig data/intermediate/C4A_complement_signature.tsv \
    --pert_info data/raw/GSE92742_Broad_LINCS_pert_info.txt \
    --out results/connectivity_by_pert.tsv
```
Aggregate across cell lines to produce drug-level reversal scores:
```bash
python src/lincs/C4A_drug_reversals.py \
    --pert-scores results/connectivity_by_pert.tsv \
    --pert-info data/raw/GSE92742_Broad_LINCS_pert_info.txt \
    --out results/C4A_drug_reversals.tsv
```

### Step 3: Train Chemprop GNN & Predict Novel Chemical Space
Format data, train the Directed Message Passing Neural Network (D-MPNN), and screen external compound libraries:
```bash
# 1. Prepare training data and generate scaffold/random split
python src/ml_chemprop/prep_and_split_for_chemprop.py \
    --input data/intermediate/C4A_training_dataset_filled.tsv \
    --out-prefix chemprop \
    --seed 42

# 2. Filter predictions (PAINS removal, Nitro exclusion, QED > 0.25, MW < 600)
python src/ml_chemprop/filter_and_cluster.py \
    --input results/chemprop_C4A_predictions.csv \
    --out results/predictions_filtered.tsv

# 3. Clean and isolate top drug-like leads
python src/ml_chemprop/clean.py \
    --input results/predictions_filtered.tsv \
    --out results/final_shortlist_chemprop.tsv
```

### Step 4: Permutation FDR Testing & CNS Feasibility
Execute a 10,000-iteration full-vector permutation null to establish empirical significance and evaluate Blood-Brain Barrier (BBB) permeability:
```bash
# Permutation null testing
python src/validation/compute_permutation_fdr.py \
    --reversals results/FINAL_OUTPUT.tsv \
    --n_perms 10000 \
    --out results/C4A_with_fdr.tsv

# CNS Feasibility (MPO and BBB scoring)
python src/validation/cns_feasibility.py \
    --input results/FINAL_OUTPUT.tsv \
    --out results/final_results_bias_adjusted.tsv
```

---

## Key Findings: Antipsychotics & Psychiatric Benchmark Hits

The unbiased transcriptomic inversion screen spontaneously recovered established psychiatric drugs, downstream signaling targets, and neuroprotective modulators (persisted in `results/antipsychotics_and_cns_candidates.tsv`):

1. **First-Generation Phenothiazine Antipsychotics:**
   - **nor-1-chlorpromazine** (Chemprop ML Probability: `0.6958`, $\text{QED} = 0.834$): Primary active metabolite and core chemical scaffold of **chlorpromazine** (the first clinical antipsychotic discovered in 1952).
   - **thioridazine** (LINCS Reversal Screen: 57 signature replicates, CNS Priority: `0.425`): Classic phenothiazine $D_2$ antagonist confirming that phenothiazines transcriptionally counter complement-mediated neuroinflammation.
2. **Dopamine Signaling Convergence Target:**
   - **GSK-3-inhibitor-IX** (Chemprop ML Probability: `0.9063`): Highest-scoring compound in the entire Chemprop model. Glycogen Synthase Kinase $3\beta$ is the primary intracellular kinase modulated downstream of dopamine $D_2$ receptor blockade and lithium action.
3. **Extrapyramidal Symptom (EPS) Antidote:**
   - **biperiden** (Chemprop ML Probability: `0.4578`, $\text{QED} = 0.832$): Central $M_1$ muscarinic antagonist specifically prescribed in clinical psychiatry to treat and reverse antipsychotic-induced extrapyramidal motor side effects.
4. **Acute Psychosis Anxiolytics & Stabilizers:**
   - **lorazepam** (`0.8441`, $\text{QED} = 0.848$), **chlordiazepoxide** (`0.8701`, $\text{QED} = 0.650$), and **alprazolam** (`0.8237`, $\text{QED} = 0.689$): Potent $\text{GABA}_\text{A}$ modulators used frontline in psychiatric emergencies to control acute agitation and reduce excitotoxic synaptic damage.
5. **Epigenetic Remodeling (Pan-HDAC Inhibitors):**
   - **panobinostat** (Reversal Score: `0.8907`), **vorinostat** (`0.2617`), and **trichostatin-a** (`0.1310`): Upregulate silenced synaptic plasticity genes (`SYN1`, `GRIN2A`) to restore cognitive integrity.
6. **Neuroendocrine & Hormone Modulators:**
   - **calcitriol** (active Vitamin $D_3$, Reversal Score: `0.1721`): Aligns with developmental vitamin D deficiency risk in schizophrenia.
   - **estrone & estradiol** (Reversal Score: `0.1427` / `0.0744`): Re-discovers the well-documented estrogen neuroprotective effect against psychotic symptoms.

Consult the full methodology report:
**[docs/METHODOLOGY_AND_ANTIPSYCHOTICS_DISCOVERY.md](docs/METHODOLOGY_AND_ANTIPSYCHOTICS_DISCOVERY.md)**

---

## Citation & References

* **Sekar, A., et al. (2016).** *Schizophrenia risk from complex variation of complement component 4.* **Nature**, 530(7589), 177–183.
* **Subramanian, A., et al. (2017).** *A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles.* **Cell**, 171(6), 1437–1452.
* **Yang, K., et al. (2019).** *Analyzing Learned Molecular Representations for Property Prediction.* **Journal of Chemical Information and Modeling**, 59(8), 3370–3388.
* **Beaulieu, J. M., et al. (2005).** *Lithium Antagonizes Dopamine-Dependent Behaviors Mediated by an AKT/Glycogen Synthase Kinase 3 Signaling Complex.* **Cell**, 122(2), 261–273.
