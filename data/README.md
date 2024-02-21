# Data
As the raw sequencing data are too large to include in this repository, this folder contains precomputed data that allows the computation of some of the figures in our manuscript.

For the Enge data we provide three files
- `engeDepthAll200V2.pkl` is a pickle file that contains a pandas dataframe with the information how many reads are aligned to each cell (i.e., the depth)
- `engeVariantsSTAR200.pkl` is a pickle file that contains a pandas dataframe with the mutation information for all cells
- `engeFilteredExpression200.h5ad` is the gene expression matrix

We also provide precalculated posterior distributions for the human data used for the figures in our manuscript.

See the iPython notebook under `./../code/geneExpressionAnalysis` for information on how to load this data.

The raw sequencing data are available from GEO or AWS and the scripts to process them with STAR are contained in the `code` subfolder.
