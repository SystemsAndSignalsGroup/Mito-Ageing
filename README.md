# Mito-Ageing
Bioinformatics pipeline for obtaining heteroplasmy features from scRNA-seq raw data. The pipeline goes from raw .fastq files and outputs standard features of a single cell workflow including aligned .bam files and expression matrices. As well as this, the heteroplasmy of mutations observed in single cells is output as a csv file which is then used for further analysis.

This code accompanies the manuscript "Cryptic mitochondrial ageing takes a lifetime and is pathophysiologically informative in single cells across tissues and species" by Alistair Green, Florian Klimm, Aidan S. Marshall, Juvid Aryaman, Patrick F. Chinnery, and Nick S. Jones.

## Prerequisites
- [STAR](https://github.com/alexdobin/STAR) for the alignment of the reads
- PySAM for variant calling
- scanpy for the gene-expression matrix analysis

## Pipeline Workflow

![Alt text](./figures/workflow/Pipeline_Workflow.png?raw=true "Title")
