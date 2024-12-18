These scripts take aligned .bam files and unfiltered expression matrices from STARsolo and performs quality control on both before provided dataframes of depth, variants and coverage for classification.
To start after running the alignment use `Expression_filtering_10x.ipynb` to perform initial quality control on the expression matrices and output a list of high quality barcodes to call variants from.
Then use `variant_submission_10x.sh` to call variants from high quality barcodes.
Finally `Mutant_classifier_10x.sh` takes the individual cells variant, depth and coverage files and performs cellular mitochondrial quality control and classification based on the number of cells any mutation appears in, as well as creates a final filtered expression matrix with all barcodes passing both quality control stages.

As these scripts need the aligned .bam files as input you cannot execute them on their own without receiving errors.
