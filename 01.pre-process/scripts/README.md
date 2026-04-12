# scripts

SLURM batch submission scripts for cluster-based preprocessing steps. These scripts run on the TSCC2 HPC and use tools such as `hicluster`, CellRanger Arc, and standard bioinformatics utilities.

## Scripts

| Script | Description |
|--------|-------------|
| `00.split_single_cell.sh` | Demultiplex Hi-C reads into per-cell contact files |
| `01.filter-contact.sh` | Filter Hi-C contacts using `hicluster` |
| `02.gene-score.sh` | Calculate gene-level accessibility scores from Hi-C data |
| `03.preproc_paired_hic_tscc2.sh` | Paired-end Hi-C read preprocessing pipeline |
| `04.proc_paired_hic_v4.sh` | Advanced Hi-C processing (alignment, dedup, contact extraction) |
| `05.compartment.sh` | Chromatin A/B compartment calling |
| `06.insulation.sh` | TAD insulation score calculation |
| `07.proc_10Xarc_RNA.sh` | CellRanger Arc RNA mapping |
| `08.proc_10Xarc_DNA.sh` | CellRanger Arc DNA (ATAC) mapping |
| `bed2bam2bw.sh` | Convert BED to BAM to BigWig format for visualization |
