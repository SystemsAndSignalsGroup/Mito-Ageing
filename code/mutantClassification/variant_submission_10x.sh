#!/bin/sh
#PBS -N call_var_sub
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-N

: '
This script will submit all the required jobs to deal with variant calling from aligned and expression matrix processed 10x data. It takes in a list of GSMs associated with the dataset as well as a list of barcodes from each GSM created by pre-filtering the expression matrix.

Make sure to edit all file paths below to match your current project

MAKE SURE TO EDIT THE NUMBER OF JOBS IN THIS SUBMISSION (-J ABOVE) TO MATCH YOUR PROJECT
'

HOME_DIR=$HOME/sc_mito_genet_pheno/Double_Alignment/Alzheimers/Metadata/ # Directory with lists of gsms and barcodes in
GSM=$(head -$PBS_ARRAY_INDEX $HOME_DIR/gsms.txt | tail -1 ) # GSM list here
BARCODE_FILE=$HOME_DIR/${GSM}_filtered_barcodes.txt # File which contains the filtered barcodes from the expression matrix
JOB_LIST=$(wc -l < $BARCODE_FILE) # Barcode list here

ERROR_FOLDER=$HOME_DIR/${GSM}_Output # Output folder here
mkdir -p $ERROR_FOLDER

# The directory of your split mitochondrial bam files (Format CB_barcode.bam)
DATA_DIR=$RDS_PROJECT/sc-ageing/ephemeral/Datasets/Alzheimers/$GSM/$GSM/MitoReads 

# Genome Name and Location
GENOME_DIR="${RDS_PROJECT}/sc-ageing/live/refs/STAR275/human"
GENOME_NAME=GRCh38.primary_assembly.genome.fa
# Genome Name (changes by species - be sure to check if using non-human data)
MITO_NAME=chrM

# Where you want your variant output files to go
OUTPUT_DIR="$HOME/sc_mito_genet_pheno/Double_Alignment/Alzheimers/Data/$GSM"
mkdir -p $OUTPUT_DIR

# Minimum depth and heteroplasmy
MIN_DEPTH=200
MIN_HF=0.05

# make sure to edit the environment called by variant_calling_10x.sh and store a local copy
SUBMISSION_SCRIPT=PATH/TO/variant_calling_10x.sh

qsub -J 1-$JOB_LIST:50 -o $ERROR_FOLDER -e $ERROR_FOLDER -V $SUBMISSION_SCRIPT
