#!/bin/sh
#PBS -N remap_test
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=16:mem=100gb
#PBS -J 1-4

set -e

module load anaconda3/personal
source activate pipeline2.0

: '
This Code will take in the STARsolo/Cellranger aligned bam file and re-align them

'
echo "$Clusters"


SAMPLE_GSM=$(head -$PBS_ARRAY_INDEX /rds/general/user/asm119/home/sc_mito_genet_pheno/Ma/brain/file_finder/unique_gsm_list.txt | tail -1 )


mkdir -p /rds/general/project/sc-ageing/ephemeral/Datasets/AidanMa/brain/$SAMPLE_GSM/$SAMPLE_GSM/remap_test

now="$(date +'%m/%d/%Y')"
now2="$(date +'%r')"
echo "$now"
echo "$now2" 

echo "Export attempt"
#export PATH=/rds/general/user/asm119/home/sc_mito_genet_pheno/samtools_install/bin:$PATH
#export PATH=/rds/general/user/asm119/home/souporcell:$PATH
export PATH=/rds/general/user/asm119/home/aidan_pipeline/aidan_pipelineV2.0/remap_package:$PATH
echo "Export success"


now="$(date +'%m/%d/%Y')"
now2="$(date +'%r')"
echo "$now"
echo "$now2" 


echo "Commencing pipeline..."
remapping.py -i /rds/general/project/sc-ageing/ephemeral/Datasets/AidanMa/brain/$SAMPLE_GSM/$SAMPLE_GSM/AlignAligned.sortedByCoord.out.bam \
-b /rds/general/project/sc-ageing/ephemeral/Datasets/AidanMa/brain/$SAMPLE_GSM/$SAMPLE_GSM/AlignSolo.out/Gene/filtered/barcodes.tsv \
-f /rds/general/project/sc-ageing/live/refs/referenceGenomeRat/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
-t 16 \
-o /rds/general/project/sc-ageing/ephemeral/Datasets/AidanMa/brain/$SAMPLE_GSM/$SAMPLE_GSM/remap_test 
echo "Completed souporcell pipeline!"


now="$(date +'%m/%d/%Y')"
now2="$(date +'%r')"
echo "$now"
echo "$now2" 

