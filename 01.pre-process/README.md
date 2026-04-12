# 01.pre-process

Preprocessing pipelines for three single-nucleus multi-omic modalities collected from human heart tissue (non-failing and failing). Each subdirectory handles one data modality, and `scripts/` contains shared SLURM submission scripts.

## Directory Structure

### [00.Multiome](./00.Multiome)
Processing of 10X Multiome (simultaneous snRNA-seq + snATAC-seq) data from two cohorts: FNIH and CAREHF.

### [01.Droplet_Paired_Tag](./01.Droplet_Paired_Tag)
Processing of Droplet Paired-Tag (DPT) data capturing histone modifications and RNA in the same nuclei.

### [02.Droplet_HiC](./02.Droplet_HiC)
Quality control and preprocessing of Droplet Hi-C data for single-cell 3D genome organization.

### [scripts](./scripts)
SLURM batch scripts for cluster-based processing steps including Hi-C contact filtering, compartment calling, CellRanger Arc mapping, and format conversion.
