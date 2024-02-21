#!/bin/sh
#PBS -N scriptHumans90
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=9:mem=600gb

set -e

module load anaconda3/personal
source activate stat-rethink

HOME_DIR=$HOME/sc_mito_genet_pheno/Analytics/monteMutInfer/
LOG_DIR=$PBS_O_WORKDIR

cp $HOME_DIR/fullMonteInfer.py .

HEALTH=human

cp $HOME_DIR/data/${HEALTH}* .

echo $NEW_INDEX
ls

echo "Starting..."
SAVER="$EPHEMERAL/Likelihoods/"
mkdir -p $SAVER
python fullMonteInfer.py -thresh 0.1 -n 90\
 -crypticDF ${HEALTH}_cryptics.pkl -cellDict ${HEALTH}_cellDict.pkl \
 -baseDict ${HEALTH}_baseDict.pkl -ageDict ${HEALTH}_ageDict.pkl \
 -logFold $LOG_DIR -saveFolder $SAVER


echo "DONE!"
