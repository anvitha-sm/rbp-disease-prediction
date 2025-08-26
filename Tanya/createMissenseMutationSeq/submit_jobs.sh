#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=02:00:00
#$ -l h_vmem=4G
#$ -N submit_mut_batches
#$ -o submitjobs_output.txt
#$ -e submitjobs_error.txt
#$ -m bea
#$ -M tchhabria@ucla.edu

INPUT_CSV="/u/home/t/tchhabri/project-kappel/gnomAD/gnomad_with_refseq_all_rows.csv"
BATCH_SIZE=5000
MAX_JOBS=500

# Count total rows
TOTAL_ROWS=$(wc -l < "$INPUT_CSV")
NUM_BATCHES=$(( (TOTAL_ROWS + BATCH_SIZE - 1) / BATCH_SIZE ))
JOB_COUNT=0

for ((i=0; i<NUM_BATCHES; i++)); do
    START=$(( i * BATCH_SIZE + 1 ))
    END=$(( (i+1) * BATCH_SIZE ))
    if [ $END -gt $TOTAL_ROWS ]; then
        END=$TOTAL_ROWS
    fi

    JOB_NAME="mut_${START}_${END}"
    OUT_FILE="mut_output_${START}_${END}.txt"
    ERR_FILE="mut_error_${START}_${END}.txt"

    # Submit individual batch
    qsub -N "$JOB_NAME" \
         -o "$OUT_FILE" \
         -e "$ERR_FILE" \
         -l h_rt=24:00:00 \
         -l h_vmem=32G \
         -m bea -M tchhabria@ucla.edu \
         /u/home/t/tchhabri/project-kappel/gnomAD/run_mutation_batch.sh "$INPUT_CSV" "$START" "$END"

    JOB_COUNT=$((JOB_COUNT + 1))
    if (( JOB_COUNT % MAX_JOBS == 0 )); then
        echo "Submitted $JOB_COUNT jobs, sleeping 10s..."
        sleep 10
    fi
done

