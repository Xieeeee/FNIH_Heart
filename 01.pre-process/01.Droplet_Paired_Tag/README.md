# 01.Droplet_Paired_Tag

Processing pipeline for Droplet Paired-Tag (DPT) data, which simultaneously profiles histone modifications and RNA from the same nuclei. This modality captures epigenetic marks (e.g., H3K27ac, H3K27me3) alongside transcriptomes.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `00.Heart_DPT_preprocess.ipynb` | Histone modification filtering and initial QC (FRiP metrics, fragment analysis) |
| `01.Heart_DPT_clustering.ipynb` | Cell clustering and UMAP visualization |
| `02.Heart_DPT_histone.ipynb` | Histone modification analysis across cell types |
| `Heart_DPT_subtype_analysis.ipynb` | Cell type subtype analysis and annotation |
