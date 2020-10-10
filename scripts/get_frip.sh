#!/bin/bash

echo -e "cell\tfrip" > qc_metrics/frip.txt

for i in */picard_bam/*.bam; do
    batch=$(echo ${i} | cut -f 1 -d/)
    cell=$(echo ${i} | rev | cut -f 1 -d/ | rev)
    total=$(samtools idxstats ${i} | addCols stdin | awk '{print $3}')
    peaks=$(intersectBed -a ${i} -b aggregate/aggregated_scATAC_peaks.narrowPeak -bed -wa -u | wc -l)
    frip=$(calc $peaks/$total | awk '{print $3*100}')
    echo -e "${batch}_${cell%_f2q30_pmd.bam}\t${frip}" >> qc_metrics/frip.txt
done
