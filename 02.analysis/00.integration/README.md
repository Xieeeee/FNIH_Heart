# 00.integration

Integration of data across multiple modalities (Multiome RNA/ATAC, Droplet Paired-Tag, Droplet Hi-C) to build a unified cell atlas of the human heart.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `00.Heart_RNA_integration.ipynb` | Integrate RNA profiles from Multiome and Droplet Paired-Tag modalities |
| `01.Heart_Hi-C_embedding.ipynb` | Generate low-dimensional embeddings from single-cell Hi-C contact matrices |
| `02.Heart_HiC_integration.ipynb` | Integrate Hi-C embeddings with RNA/ATAC-based cell annotations |
| `03.plot_integration.ipynb` | Visualize integration results (UMAPs, cross-modality comparisons) |
