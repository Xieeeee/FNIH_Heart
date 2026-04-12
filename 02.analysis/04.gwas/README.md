# 04.gwas

GWAS enrichment and fine-mapping analyses to map genetic associations onto cell type-specific regulatory programs. Combines LD score regression (LDSC) for heritability partitioning with direct variant-to-gene mapping through Hi-C loops.

## Subdirectories

### LDSC/
LD Score Regression for heritability enrichment in cell type-specific cCREs.

| Notebook | Description |
|----------|-------------|
| `1_Munge_GWAS.ipynb` | Format and munge GWAS summary statistics for LDSC input |
| `2_RunLDSC_V2.ipynb` | Run partitioned LDSC heritability analysis |
| `3_PlotResults_V5.ipynb` | Visualize LDSC enrichment results across cell types and traits |

### LDSC_Chromatin_states/
LDSC analysis stratified by chromatin states (e.g., active enhancer, promoter, repressed).

| Notebook | Description |
|----------|-------------|
| `2_RunLDSC_V2.ipynb` | Run LDSC partitioned by chromatin state annotations |
| `3_PlotResults_V5.ipynb` | Visualize chromatin state-specific heritability enrichment |

### Traits_intersect/
Fine-mapping and variant-to-gene linking using Hi-C loops and GRN data.

| Notebook | Description |
|----------|-------------|
| `1_MakeBedTraits-CHROMBP.ipynb` | Create BED files for GWAS trait-associated loci |
| `1_Traits_PeaksCalls_intersect_filtered_SNPs-V2.ipynb` | Intersect GWAS fine-mapped variants with cell type-specific cCREs |
| `GSEA_FORA_GRN.ipynb` | Gene set enrichment and over-representation analysis on GRN target genes |
| `HF_FineMap_analysis_hg38_4chambersDiff_HicLoops_allSNPs_AndFilterSNPs--V8.ipynb` | Fine-mapping of heart failure variants using Hi-C loops to link SNPs to target genes across chambers |
