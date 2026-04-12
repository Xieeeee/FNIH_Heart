# scripts

Shared R utility functions used across analysis notebooks.

## Files

| Script | Description |
|--------|-------------|
| `basics.R` | Core utilities: color palettes (kanagawa, ukiyo, etc.), publication ggplot2 theme, Seurat RNA processing wrappers (`RunRNA`, `RunXPM`), DNA bin masking (`RunDNAMask`), barcode filtering, Jaccard similarity clustering, and sparse matrix operations |
| `DPT_help.R` | Droplet Paired-Tag helpers: mouse-to-human gene conversion, FRiP calculation and visualization (`ImportArcFRiP`, `PlotArcFRiP`), fragment length distribution plots, and multiome barcode pairing (`PairArc`) |
| `snapATAC2.JSStest.R` | Command-line tool for Jensen-Shannon divergence-based differential peak calling from snapATAC2 h5ad objects with empirical FDR estimation |
