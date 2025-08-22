#!/bin/bash

# Count rows (excluding header)
#total_rows=$(tail -n +2 Combined_Chunks.csv | wc -l)
total_rows=(130)
#echo "Total indices available: $total_rows"
# Submit one job per row
for row in "${total_rows[@]}"; do
    echo "Resubmitting job for row $row"
    qsub test.sh "$row"
done

