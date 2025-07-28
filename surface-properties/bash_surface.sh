#!/bin/bash
#$ -cwd
#$ -j y
#$ -o pymol_analysis_log.$JOB_ID
#$ -e error.$JOB_ID
#$ -l h_rt=6:00:00,h_data=3G
#$ -pe shared 12
#$ -t 1-35
#$ -M $USER@mail
#$ -m be

echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started at: $(date)"
echo ""

. /u/local/Modules/default/init/modules.sh
module load miniforge
conda activate myconda
python surface_analysis.py

echo ""
echo "Job $JOB_ID ended on: $(hostname -s)"
echo "Job $JOB_ID ended at: $(date)"

