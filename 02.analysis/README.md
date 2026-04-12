# 02.analysis

Downstream analysis workflows and code for reproducing manuscript figures. Analyses span multi-modal data integration, regulatory element identification, gene regulatory network inference, and genetic association mapping.

## Directory Structure

### [00.integration](./00.integration)
Cross-modality data integration (RNA, ATAC, Hi-C) using Harmony and custom embedding approaches.

### [01.cCRE](./01.cCRE)
Candidate cis-regulatory element (cCRE) analysis, including differential accessibility (ATAC) and expression (RNA) between heart failure and control using DESeq2.

### [02.GRN](./02.GRN)
Gene regulatory network (GRN) inference linking transcription factors to target genes via cCRE connections.

### [03.Hi-C](./03.Hi-C)
Compartment and chromatin domain analysis for Droplet Hi-C data

### [04.gwas](./04.gwas)
GWAS enrichment analysis using LD score regression (LDSC) and fine-mapping of causal variants through integration with cell type-specific regulatory maps and Hi-C loops.

### [scripts](./scripts)
Shared R utility functions.
