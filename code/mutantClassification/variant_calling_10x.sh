#!/bin/sh
#PBS -N call_variant
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb

: '
The only thing that should need to be edited in this file is the name of your anaconda environment 
'

set -e
module load anaconda3/personal
source activate align # SET YOUR NAME HERE


: '
No need to edit below here
'

#address to the Ali's variant calling python script
THE_PYTHON_SCRIPT="$RDS_PROJECT/sc-ageing/live/ali_pipeline/"

for (( i=$PBS_ARRAY_INDEX; i<$PBS_ARRAY_INDEX+50; i++ ))
do
#read in the cell file barcodes and index the array jobs by each cell
CELL=$(head -$i  $BARCODE_FILE | tail -1 )
CELL_FILE=CB_${CELL}.bam
echo $CELL
echo $CELL_FILE

if [[ -e $DATA_DIR/$CELL_FILE ]]; then

#exit immediatley if the file fails
cp $DATA_DIR/${CELL_FILE} .
samtools index ${CELL_FILE}

cp $THE_PYTHON_SCRIPT/snp_del_ins_caller_coverag.py .

python snp_del_ins_caller_coverage.py -bam_file_name ${CELL_FILE} \
  -mito_file $GENOME_DIR/$GENOME_NAME -sample_id ${GSM}_$CELL -min_depth $MIN_DEPTH \
  -min_hf $MIN_HF -mito_name $MITO_NAME

echo "variant calling done"
[[ -e ${GSM}_${CELL}_variants.csv ]] && cp ${GSM}_${CELL}_variants.csv $OUTPUT_DIR
cp ${GSM}_${CELL}_depth.csv $OUTPUT_DIR
cp ${GSM}_${CELL}_coverage.pkl $OUTPUT_DIR
echo "Copied Across!"
else
echo "$CELL_FILE not found - moving on"
fi

done
echo "Finished!"
