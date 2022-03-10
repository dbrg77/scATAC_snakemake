for i in */picard_bam/*.bam; do
    batch=$(echo ${i} | cut -f 1 -d/)
    fn=$(echo ${i} | rev | cut -f 1 -d/ | rev)
    cell=${fn%_f2q30_pmd.bam}
    samtools sort -@ 12 -n -T ${cell}_tmp ${i} | \
    bamToBed -bedpe -i - | \
    awk -v cn="${batch}_${cell}" 'BEGIN{OFS="\t"}{ if($9=="+"){print $1, $2, $6, cn, "1"} else {print $1, $5, $3, cn, "1"} }'
done
