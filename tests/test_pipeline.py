import os
import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestPipelineImports(unittest.TestCase):
    """Ensure all core pipeline modules can be imported without syntax or import errors."""

    def test_imports(self):
        import src.signature.build_c_4_a_complement_signature
        import src.signature.convert_signature_to_entr
        import src.lincs.analyse_reversal
        import src.lincs.compute_connectivity_lincs
        import src.lincs.C4A_drug_reversals
        import src.lincs.extract_lincs_metadata
        import src.lincs.gctx
        import src.lincs.map_to_lincs
        import src.ml_chemprop.clean
        import src.ml_chemprop.filter_and_cluster
        import src.ml_chemprop.prep_and_split_for_chemprop
        import src.ml_chemprop.reversal
        import src.validation.cns_feasibility
        import src.validation.compute_permutation_fdr
        import src.validation.compute_reversal_significance
        import src.validation.finalize_shortlist
        import src.validation.mapping_audit
        import src.validation.prioritisation
        import src.validation.summarize_signature_consistency
        import src.utilities.check_genes
        import src.utilities.fetch_smiles_pubchem
        import src.utilities.map_genes
        import src.utilities.map_lincs_to_smiles
        import src.utilities.quick_sanity


class TestSignatureLogic(unittest.TestCase):
    """Test column detection and combined z-score logic."""

    def test_detect_columns(self):
        from src.signature.build_c_4_a_complement_signature import detect_columns
        df = pd.DataFrame({
            "gene_symbol": ["C4A", "C3", "C2"],
            "logFC": [1.5, -0.8, 0.4],
            "P.Value": [0.01, 0.05, 0.2]
        })
        gene_col, logfc_col, p_col = detect_columns(df)
        self.assertEqual(gene_col, "gene_symbol")
        self.assertEqual(logfc_col, "logFC")
        self.assertEqual(p_col, "P.Value")


class TestFDRAndStatistics(unittest.TestCase):
    """Test Benjamini-Hochberg and one-sided permutation test statistics."""

    def test_benjamini_hochberg(self):
        from src.validation.compute_permutation_fdr import benjamini_hochberg
        pvals = [0.001, 0.01, 0.04, 0.05, 0.5]
        qvals = benjamini_hochberg(pvals)
        self.assertEqual(len(qvals), len(pvals))
        # Monotonicity check
        for i in range(len(qvals) - 1):
            self.assertLessEqual(qvals[i], qvals[i + 1] + 1e-9)

    def test_one_sided_p_values_from_null(self):
        from src.validation.compute_permutation_fdr import one_sided_p_values_from_null
        obs = [-2.5, 0.5]
        # Null matrix shape: (1000 perms, 2 observations)
        np.random.seed(42)
        null = np.random.normal(0, 1, size=(1000, 2))
        pvals = one_sided_p_values_from_null(obs, null, direction="negative")
        self.assertLess(pvals[0], 0.05)  # Extreme negative should be significant
        self.assertEqual(pvals[1], 1.0)  # Positive observed has p=1 in negative direction


class TestResultsIntegrity(unittest.TestCase):
    """Verify that core output datasets exist and contain required biological data."""

    def test_benchmark_antipsychotics_present(self):
        tsv_path = PROJECT_ROOT / "results" / "antipsychotics_and_cns_candidates.tsv"
        self.assertTrue(tsv_path.exists(), f"Missing {tsv_path}")
        df = pd.read_csv(tsv_path, sep="\t")
        self.assertGreater(len(df), 5)
        self.assertIn("Compound", df.columns)

        compounds = df["Compound"].str.lower().tolist()
        # Verify key benchmark rediscoveries
        self.assertTrue(any("chlorpromazine" in c for c in compounds))
        self.assertTrue(any("gsk-3" in c or "gsk3" in c for c in compounds))

    def test_final_output_exists(self):
        tsv_path = PROJECT_ROOT / "results" / "FINAL_OUTPUT.tsv"
        self.assertTrue(tsv_path.exists(), f"Missing {tsv_path}")
        df = pd.read_csv(tsv_path, sep="\t")
        self.assertGreater(len(df), 800)
        self.assertIn("reversal_score", df.columns)


if __name__ == "__main__":
    unittest.main()
