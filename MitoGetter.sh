#!/bin/bash

: '
This Script will take in a full bam file which has been been tagged with CB, the cell barcodes
It will then create a BAM file per cell, containing only the mitochondrial reads
'


GENOME_DIR="${RDS_PROJECT}/sc-ageing/live/refs/referenceGenomeRat"
#GENOME_NAME=GRCh38.primary_assembly.genome.fa


#The arguments the script needs
#Get the .bam file to be worked on
BAMFILE=$1
#Get the name of the mitochondrial chromosome
ChrMt=$2
#Output Directory
Output=$3
#barcode whitelist used
whitelist=$4

: ' STEP 1 - Getting only the mitochondrial reads
'
samtools index $BAMFILE
samtools view -b $BAMFILE $ChrMt > ${ChrMt}.bam
samtools sort -t CB -o ${ChrMt}_sorted_tags.bam ${ChrMt}.bam

# copy the files back to the output directory
cp ${ChrMt}_sorted_tags.bam $Output/
cp ${ChrMt}.bam $Output/${ChrMt}.bam
echo "	Mitoreads obtained"


# command to keep only the reads barcodes
samtools view -D CB:$whitelist -o sorted_tags_MT_only_BC.bam ${ChrMt}_sorted_tags.bam
cp sorted_tags_MT_only_BC.bam $Output

echo "	splitting bam per cell"
: ' Splitting the BAM file

It calls a python script called split_script.py to split the BAM according to the aligned cell BARCODE tags. This file only requires the location of the BAM file to be split as an argument
'

mkdir -p MitoReads
split_script.py sorted_tags_MT_only_BC.bam MitoReads
cp -r MitoReads $Output/MitoReads

#make a list of individual cell file names
ls MitoReads > cells.txt
#copy this list out for later viewing
cp cells.txt $Output
echo "	Split success"
