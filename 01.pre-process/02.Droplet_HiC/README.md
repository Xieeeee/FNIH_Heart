# 02.Droplet_HiC

Preprocessing of Droplet Hi-C data for single-cell 3D genome organization profiling. Hi-C contact matrices are used to derive chromatin compartments, TAD boundaries, and loop structures.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `00.DHC_QC_preprocess.ipynb` | Quality control and preprocessing of Droplet Hi-C libraries |

Additional Hi-C processing steps (contact filtering, compartment calling, insulation scoring) are handled by the SLURM scripts in [`../scripts/`](../scripts/).
