# Mito-Ageing


Bioinformatics pipeline for obtaining heteroplasmy features from scRNA-seq raw data. The pipeline goes from raw .fastq files and outputs standard features of a single cell workflow including aligned .bam files and expression matrices. As well as this, the heteroplasmy of mutations observed in single cells is output as a csv file which is then used for further analysis.

This code accompanies the manuscript "Cryptic mitochondrial ageing takes a lifetime and is pathophysiologically informative in single cells across tissues and species" by Alistair Green, Florian Klimm, Aidan S. Marshall, Juvid Aryaman, Patrick F. Chinnery, and Nick S. Jones.

## Prerequisites
- [STAR](https://github.com/alexdobin/STAR) for the alignment of the reads
- PySAM for variant calling
- scanpy for the gene-expression matrix analysis

## Pipeline Workflow

### How-to

This project consist of four steps
1. Processing of raw sequencing reads
- Downloading the data from GEO or AWS
- Alignment of reads to reference genome with STAR
- Filtering out mitochondrial reads for variant calling
2. Characterisation of mitochondrial mutations at a single-cell level
- Quality control of mitochondrial mutation data
- Comparing mutations across the data set to identify `cryptic` mutations
3. Figure creating for analysis of mitochondrial mutations
- Comparison with pseudo-bulk heteroplasmy
- Computation of site-frequency spectrum
- Selection effects
- Increasing homoplasmies with age
- Increasing difference between SFS with increasing age difference of donors
4. Gene-expression analysis
- Quality control of gene-expresison matrix
- Incorporation of mitochondrial


Due to space limitations, we do not provide raw sequencing read data in this GitHub (but they can be downloaded from GEO and AWS). Rather, we provide preprocessed mutation data and gene-expression matrix for one example data for the human pancreas [link](https://pubmed.ncbi.nlm.nih.gov/28965763/) .


![Alt text](./figures/workflow/Pipeline_Workflow.png?raw=true "Title" =250x)
