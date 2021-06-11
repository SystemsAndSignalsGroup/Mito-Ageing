We provide example scripts for different steps in our analysis pipeline.

1. Processing of raw sequencing reads: `./alignment`
2. Characterisation of mitochondrial mutations at a single-cell level: `./mutantClassification`
3. Figure creating for site-frequency-spectrum and other statistics: `./variantAnalysis`
4. Gene-expression analysis: `./geneExpressionAnalysis`

For steps (3) and (4) we provide precomputed data such that figures from our manuscript can be reproduced. For step (1) you first have to download the raw sequencing reads from GEO or AWS, respectively. The output from step (1) can then be used as input for step (2), which produces data that can be used for further analysis.
