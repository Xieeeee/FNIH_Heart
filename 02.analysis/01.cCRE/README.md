# 01.cCRE

Candidate cis-regulatory element (cCRE) analysis and differential testing between heart failure (HF) and non-failing controls (CTR). This directory identifies cell type-specific regulatory elements and quantifies disease-associated changes in both chromatin accessibility and gene expression.

## Notebooks

- **`cCRE_analysis.ipynb`** - Main cCRE analysis: identification, annotation, and characterization of candidate cis-regulatory elements across cell types and chambers.

## Subdirectories

### DESEQ_ATAC/
Differential chromatin accessibility analysis using DESeq2.

- **`1_Dowsntream_files/`**
  - `0_Create_CellTypeSpecificPeaks.ipynb` - Identify cell type-specific ATAC peak sets
  - `1_4Chambers_ATAC_Dowstream_Files-donor.ipynb` - Create per-donor ATAC count matrices across all four chambers
- **`2_Deseq/`**
  - `3.1_4Chambers_DEseq2_perdonor_HFvsCTR.ipynb` - DESeq2 differential accessibility (HF vs. CTR) per donor, per chamber

### DESEQ_RNA/
Differential gene expression analysis using DESeq2.

- **`1_Dowsntream_files/`**
  - `1_4Chambers_RNA_Dowstream_Files-donor.ipynb` - Create per-donor RNA pseudobulk count matrices across all four chambers
- **`2_Deseq/`**
  - `3.1_4Chambers_DEseq2_perdonor_HFvsCTR.ipynb` - DESeq2 differential expression (HF vs. CTR) per donor, per chamber
