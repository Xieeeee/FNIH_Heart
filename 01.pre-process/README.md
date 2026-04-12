# 01.pre-process

Preprocessing pipelines for three single-nucleus multi-omic modalities collected from human heart tissue (non-failing and failing). Each subdirectory handles one data modality, and `scripts/` contains shared SLURM submission scripts.

## Directory Structure

### [00.Multiome](./00.Multiome)
Processing of 10X Multiome data (snATAC + RNA) from multiplexed and single-channel reactions.

### [01.Droplet_Paired_Tag](./01.Droplet_Paired_Tag)
Processing of Droplet Paired-Tag (DPT) data capturing histone modifications and RNA in the same nuclei.

### [02.Droplet_HiC](./02.Droplet_HiC)
Quality control and preprocessing of Droplet Hi-C data for single-cell 3D genome organization.

### [scripts](./scripts)
SLURM batch scripts for Droplet Hi-C and pre-processing.
