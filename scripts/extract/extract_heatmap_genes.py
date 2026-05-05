"""
Extract EdgeR logFC and FDR for heatmap genes across selected cell types.

Input:  data/EdgeR/OIL_vs_CORT/<cell_type>_edgeR_results.tsv
Output: data/processed/heatmap_neuron_glia.csv
"""

from pathlib import Path
import pandas as pd

EDGER_DIR = Path("data/EdgeR/OIL_vs_CORT")
OUT_CSV   = Path("data/processed/heatmap_neuron_glia.csv")

# ── Gene manifest by theme ─────────────────────────────────────────────────────
THEMES = {
    "1. GR desensitization":   ["Nr3c1", "Nr3c2", "Sgk1", "Fkbp5", "Ddit4"],
    "2. Neurotrophic collapse": ["Bdnf", "Ntrk2", "Ptn", "Sort1"],
    "3. Circadian disruption":  ["Per1", "Per2", "Per3", "Cry1"],
    "4. Neuroinflammation":     ["Il6r", "Nfkbia", "Tsc22d3", "Ifngr1", "Mif"],
    "5. E/I balance":           ["Uchl1", "Fstl5", "Chl1"],
    "7. cAMP/CREB plasticity":  ["Adcy8", "Nr4a3", "Creb5", "Pde9a"],
    "8. Wnt/synaptic adhesion": ["Wls", "Pcdh15", "Pcdh7", "Uqcrq"],
}

# ── Cell type manifest (label → filename, class) ───────────────────────────────
CELL_TYPES = {
    "CA1-ProS":  ("016_CA1-ProS_Glut_edgeR_results.tsv",  "Hippocampal neuron"),
    "CA3":       ("017_CA3_Glut_edgeR_results.tsv",        "Hippocampal neuron"),
    "DG":        ("037_DG_Glut_edgeR_results.tsv",         "Hippocampal neuron"),
    "L2/3 IT":   ("007_L2_3_IT_CTX_Glut_edgeR_results.tsv",  "Cortical glutamatergic"),
    "L4/5 IT":   ("006_L4_5_IT_CTX_Glut_edgeR_results.tsv",  "Cortical glutamatergic"),
    "L5 ET":     ("022_L5_ET_CTX_Glut_edgeR_results.tsv",    "Cortical glutamatergic"),
    "L6 IT":     ("004_L6_IT_CTX_Glut_edgeR_results.tsv",    "Cortical glutamatergic"),
    "L6 CT":     ("030_L6_CT_CTX_Glut_edgeR_results.tsv",    "Cortical glutamatergic"),
    "Pvalb":     ("052_Pvalb_Gaba_edgeR_results.tsv",      "Interneuron"),
    "Sst":       ("053_Sst_Gaba_edgeR_results.tsv",        "Interneuron"),
    "Astro-NT":  ("318_Astro-NT_NN_edgeR_results.tsv",     "Non-neuronal"),
    "Astro-TE":  ("319_Astro-TE_NN_edgeR_results.tsv",     "Non-neuronal"),
    "OPC":       ("326_OPC_NN_edgeR_results.tsv",          "Non-neuronal"),
    "Oligo":     ("327_Oligo_NN_edgeR_results.tsv",        "Non-neuronal"),
    "Microglia": ("334_Microglia_NN_edgeR_results.tsv",    "Non-neuronal"),
    "Endo":      ("333_Endo_NN_edgeR_results.tsv",         "Non-neuronal"),
    "Peri":      ("331_Peri_NN_edgeR_results.tsv",         "Non-neuronal"),
    "VLMC":      ("330_VLMC_NN_edgeR_results.tsv",         "Non-neuronal"),
    "CHOR":      ("325_CHOR_NN_edgeR_results.tsv",         "Non-neuronal"),
}

rows = []
for cell_type, (filename, cell_class) in CELL_TYPES.items():
    df = pd.read_csv(EDGER_DIR / filename, sep="\t", index_col=0)
    for theme, genes in THEMES.items():
        for gene in genes:
            if gene in df.index:
                r = df.loc[gene]
                rows.append({
                    "gene":       gene,
                    "theme":      theme,
                    "cell_type":  cell_type,
                    "cell_class": cell_class,
                    "logFC":      r["logFC"],
                    "FDR":        r["FDR"],
                })
            else:
                rows.append({
                    "gene":       gene,
                    "theme":      theme,
                    "cell_type":  cell_type,
                    "cell_class": cell_class,
                    "logFC":      float("nan"),
                    "FDR":        float("nan"),
                })

out = pd.DataFrame(rows)
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_CSV, index=False)
print(f"Saved {OUT_CSV}  ({len(out)} rows, {out['gene'].nunique()} genes, {out['cell_type'].nunique()} cell types)")