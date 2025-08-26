#!/bin/bash
#$ -cwd
#$ -V
#$ -m bea
#$ -M tchhabria@ucla.edu

# Runtime and memory for individual batch
#$ -l h_rt=24:00:00
#$ -l h_vmem=32G

# Arguments
INPUT_CSV="$1"
START="$2"
END="$3"

echo "Running batch $START to $END"
date

# Load conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pymol-env

# Run Python script for this batch
python /u/home/t/tchhabri/project-kappel/gnomAD/run_mutations.py \
    --input "$INPUT_CSV" \
    --start "$START" \
    --end "$END"

date
echo "Finished batch $START to $END"
