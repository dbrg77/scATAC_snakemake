for i in */picard_bam/*.bam; do
    fn=$(echo ${i} | rev | cut -f 1 -d/ | rev)
    cell=${fn%_f2q30_pmd.bam}
    samtools view -f 35 ${i} | awk -v cn=${cell} 'BEGIN{OFS="\t"}{print $3, $4-1, $4+$9, cn, "1"}'
done
