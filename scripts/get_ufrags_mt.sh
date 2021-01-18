#!/bin/bash

echo -e "cell\tuniq_nuc_frags" > qc_metrics/uniq_nuc_frags.txt
echo -e "cell\tmt_content" > qc_metrics/mt_content.txt

for i in */picard_bam/*.bam; do
    batch=$(echo ${i} | cut -f 1 -d/)
    cell=$(echo ${i} | rev | cut -f 1 -d/ | rev)
    tread=$(samtools idxstats ${i} | addCols stdin | awk '{print $3}')
    mt=$(samtools idxstats ${i} | grep chrM | awk '{print $3}')
    nread=$(calc ${tread}-${mt} | awk '{print $3}')
    nfrag=$(calc ${nread}/2 | awk '{print $3}' | cut -f 1 -d '.')
    p=$(calc ${mt}/${tread}*100 | awk '{print $3}')
    echo -e "${batch}_${cell%_f2q30_pmd.bam}\t${nfrag}" >> qc_metrics/uniq_nuc_frags.txt
    echo -e "${batch}_${cell%_f2q30_pmd.bam}\t${p}" >> qc_metrics/mt_content.txt
done
