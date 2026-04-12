# 00.Multiome

Processing pipeline for 10X Multiome data (simultaneous snRNA-seq + snATAC-seq from the same nuclei). Two independent cohorts are processed separately before downstream integration.

## Notebooks

- **`01.CAREHF_clustering.ipynb`** - Clustering and QC for the CAREHF Multiome cohort (single-channel processed data).

## Subdirectories

### Processing_CAREHF/
Step-by-step preprocessing for the CAREHF cohort:

| Notebook | Description |
|----------|-------------|
| `1_PrefilterQC.ipynb` | Initial quality control and filtering |
| `2_Switch_to_ReferenceATACpeaks.ipynb` | Re-quantify ATAC using a unified reference peak set |
| `3_SoupX.ipynb` | Ambient RNA removal with SoupX |
| `4.1_Amulet.ipynb` | Doublet detection using AMULET (DNA-based) |
| `4.2_Scrublet.SamplePrep.ipynb` | Prepare inputs for Scrublet doublet detection |
| `4.3_Scrublet_Run.ipynb` | Run Scrublet doublet detection |
| `4.4_doublet.cleanup-either.ipynb` | Remove doublets flagged by either AMULET or Scrublet |
| `5_CAREHF_Upload.ipynb` | Finalize and prepare data for upload |

### Processing_MultiomeFNIH/
Full preprocessing pipeline for the main FNIH Multiome cohort:

| Notebook | Description |
|----------|-------------|
| `1_Merge_Filter_rmvMT_Harmonize[lane]_annotate.ipynb` | Merge lanes, filter cells, remove mitochondrial reads, Harmony batch correction, initial annotation |
| `2_Sil_QC-V2-Compartment.ipynb` | Silhouette-based QC and chromatin compartment analysis |
| `3.2_CleanupMap.ipynb` | Additional cleanup and re-mapping |
| `4_Finalmap_annotation.ipynb` | Final cell type annotation |
| `5_Finalize_objects.ipynb` | Create finalized Seurat objects |
| `5.1_CallPeaks_ByChamber_SubTypes.ipynb` | Call ATAC peaks stratified by chamber and cell subtype |
| `6_Papers_Plots.ipynb` | Generate publication figures |
| `7_FNIH_Multiome_FinalizeForUpload.ipynb` | Prepare final objects for public data upload |
| `8_CreateAnnData_Object.ipynb` | Convert Seurat objects to Python AnnData format |
