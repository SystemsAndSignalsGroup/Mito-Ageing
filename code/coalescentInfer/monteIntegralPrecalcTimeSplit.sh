#!/bin/sh
#PBS -N monteMake
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -J 0-799

set -e

cp monte_integral_timeArray.py .
echo "Starting..."

FID=$PBS_ARRAY_INDEX

echo $FID

NUMBA_NUM_THREADS=8 python monte_integral_timeArray.py -threads 8 -time_index $FID -tot_levels 90 \
        -save_path ./Likelihoods

echo "DONE!"
