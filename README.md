# Visualisation for NGS Data


## Installation guide

The scripts can be run directly from source code.

## Per base coverage plot

Creates a coverage plot for each nucleotide. Specifically useful for viruses. Can read in multi fasta data as a [reference genome](https://github.com/jonas-fuchs/General-NGS-visualisation-tools/tree/main/Coverage%20plot/Reference.fasta) and generates mfacet plots based on Qualimaps per base [coverage output](https://github.com/jonas-fuchs/General-NGS-visualisation-tools/tree/main/Coverage%20plot/Coverage.tsv).

```R
Coverage plot/Coverage.R
```

![Example output](https://github.com/jonas-fuchs/General-NGS-visualisation-tools/tree/main/Coverage%20plot/Coverage.png)