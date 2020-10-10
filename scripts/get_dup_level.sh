#!/bin/bash

echo -e "cell\tdup_level" > qc_metrics/dup_level.txt

for i in */picard_log/*.out; do
    dup=$(grep Unknown ${i} | cut -f 9)
    batch=$(echo ${i} | cut -f 1 -d/)
    cell=$(echo ${i} | rev | cut -f 1 -d/ | rev)
    echo -e "${batch}_${cell%_f2q30_pmd.out}\t${dup}" >> qc_metrics/dup_level.txt
done
