#!/bin/bash

# === CONFIGURATION ===
INPUT_DIR="/u/project/kappel/fraza/CIF Files/Biological_Assemblies"
BATCH_SIZE=25
SCRIPT_DIR="/u/home/a/anvitha/project-kappel"
QSUB_SCRIPT="$SCRIPT_DIR/bash_surface.sh"

# === PREP ===
mkdir -p batch_lists
mkdir -p logs

# === LIST FILES SAFELY ===
# Write full list of filenames (basename only) to a master list
find "$INPUT_DIR" -maxdepth 1 -type f -name '*.cif' -exec basename {} \; | sort > all_files.txt

TOTAL=$(wc -l < all_files.txt)
NUM_BATCHES=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Total files: $TOTAL"
echo "Batch size: $BATCH_SIZE"
echo "Total batches: $NUM_BATCHES"

# === SPLIT INTO BATCHES ===
split -l $BATCH_SIZE all_files.txt batch_lists/batch_

# === SUBMIT JOBS ===
for file in batch_lists/batch_*; do
    # Rename for clarity
    batch_id=$(basename "$file" | cut -d_ -f2)
    new_file="batch_lists/batch_${batch_id}.txt"
    mv "$file" "$new_file"

    echo "Submitting $new_file with $(wc -l < "$new_file") files"
    qsub "$QSUB_SCRIPT" "$new_file"
done

