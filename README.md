# kTYPr_EcoEpidem

This repository contains code and data used for the analysis of *E. coli* transporter-dependent-capsule (K-type) diversity, ecology and epidemiology. 

## Profiled collections:

```text
- RefSeq inputdir/refseq37_flank_20260107_ktyps.tsv
- NCBI inputdir/ncbi_flank_20260108_ktyps.tsv
- Egli inputdir/egli_flank_20260115_ktyps.tsv
- Gut MAGs inputdir/commensal_flank_20260115_ktyps.tsv
```

## Main repository structure

This refers to Figures 3-4, Extended Data Figures 3, 6-9 and Tables 12-14  in the study "The natural diversity of *E. coli* transporter-dependent capsules".

```text
code_Ecocapsules.Rmd   # Main analysis and figure-generation code
inputdir/              # Input data
outputdir/             # Processed outputs including Extended Data Tables
plotdir/	       # Empty directory to reproduce all figures.
```

## Additional analysis folder structure

The [additional_analysis](additional_analysis) folder contains documentation and scripts to run the presented ANI dereplication analyses, BLASTN/BLASTP (Extended Data Table 6), kTYPr-Kaptive comparative and perturbation analysis (Extended Data Table 17).
