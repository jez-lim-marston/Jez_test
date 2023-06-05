# Snakemake workflow

## Title: Comparing HERV Expression Between Responder and Non-responder to Immunotherapy in Metastatic Melanoma

Bulk RNA-Seq Telescope analysis of 2 metastatic melanoma datasets. Input the metadata table that include the accessions needed to download the raw sequencing files, output gene/HERV count matrices.

## Dataset Citations

1. Hugo, W., Zaretsky, J. M., Sun, L. U., Song, C., Moreno, B. H., Hu-Lieskovan, S., ... & Lo, R. S. (2016). Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell*, 165(1), 35-44.

2. Van Allen, E. M., Miao, D., Schilling, B., Shukla, S. A., Blank, C., Zimmer, L., ... & Garraway, L. A. (2015). Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. *Science*, 350(6257), 207-211.

## Workflow Graphs

### To get DAG:
 
```snakemake --profile profiles/aws  --forceall --dag | dot -Tpdf > dag.pdf```

### To get rule graph:

```snakemake --profile profiles/aws  --forceall --rulegraph | dot -Tpdf > rulegraph.pdf```

### To get file graph:

```snakemake --profile profiles/aws  --forceall --filegraph | dot -Tpdf > filegraph.pdf```

### To run pipeline:

```snakemake --profile profiles/aws/ all```

## To modify pipeline:

Change sample download table and method. This pipeline uses ```prefetch```then```fastq-dump``` to download files.


## Usage

 If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository: https://github.com/nixonlab/Mets_CM_bulk
