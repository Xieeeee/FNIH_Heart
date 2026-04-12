# 03.Hi-C

Hi-C chromatin conformation analysis for the FNIH Heart dataset. This directory covers visualization of contact matrices, A/B compartment analysis, topological domain boundary calling, and differential Hi-C compartment (DHC) analysis between heart failure (HF) and non-failing controls.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `00.view_hic.ipynb` | Visualize Hi-C contact matrices at multiple resolutions across cell types (vCM, Fibroblast, Endothelial, Myeloid, Lymphoid) and HF conditions; includes imputed and raw maps, rotated triangle views, and HF-vs-control subtraction maps |
| `01.compartment_analysis.ipynb` | A/B compartment calling using bulk PCA models, cell type-specific compartment scores, saddle plots, compartment strength quantification, differential compartment analysis with dcHiC, contact decay curves by compartment, and single-cell compartment scores |
| `02.domain_analysis.ipynb` | Topological domain boundary calling at donor-pseudobulk level, boundary saturation analysis, differential boundary detection across cell types using chi-squared tests, insulation score computation, and boundary pileup analysis comparing HF vs control |
| `04.DHC_analysis` | Differential Hi-C compartment analysis in R: defines `diffGAD` function for Wilcoxon-based differential compartment testing across cell types and conditions |
