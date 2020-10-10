configfile: 'config.json'

from glob import iglob
import pandas as pd

r1 = iglob('*/fastq/*_r1.fq.gz')

samples = pd.DataFrame()
samples['r1'] = [i for i in r1]
samples[['rep', 'cell']] = samples.r1.str.extract('(.*)/fastq/(.*)_r1\.fq\.gz', expand=True)
replicate = list(set(samples['rep'].values))

rule all:
    input:
        list(samples['rep'] + '/fastp_fq/' + samples['cell'] + '_r1_trimmed.fq.gz'),
        list(samples['rep'] + '/fastp_fq/' + samples['cell'] + '_r2_trimmed.fq.gz'),
        list(samples['rep'] + '/fastp_report/' + samples['cell'] + '_fastp.html'),
        list(samples['rep'] + '/fastp_report/' + samples['cell'] + '_fastp.json'),
        list(samples['rep'] + '/hisat2_mapped/' + samples['cell'] + '_f2q30.bam'),
        list(samples['rep'] + '/hisat2_log/' + samples['cell'] + '_aln_sum.txt'),
        list(samples['rep'] + '/picard_bam/' + samples['cell'] + '_f2q30_pmd.bam'),
        list(samples['rep'] + '/picard_bam/' + samples['cell'] + '_f2q30_pmd.bam.bai'),
        list(samples['rep'] + '/picard_log/' + samples['cell'] + '_f2q30_pmd.out'),
        list(samples['rep'] + '/isize_hist/' + samples['cell'] + '_isize.hist'),
        ['{}/bam_file_list.txt'.format(i) for i in replicate],
        ['{}/f2q30_merged.bam'.format(i) for i in replicate],
        'aggregate/f2q30_merged.bam',
        'aggregate/f2q30_merged_pmd.bam',
        'aggregate/f2q30_merged_pmd.bed',
        'aggregate/f2q30_merged_pmd.out',
        'aggregate/f2q30_merged_pmd_isize.hist',
        'aggregate/aggregated_scATAC_peaks.xls',
        'aggregate/aggregated_scATAC_peaks.narrowPeak',
        'aggregate/aggregated_scATAC_treat_pileup.bdg',
        'aggregate/aggregated_scATAC_control_lambda.bdg',
        'aggregate/aggregated_scATAC_treat_pileup.bw',
        list(samples['rep'] + '/count/' + samples["cell"] + '.count'),
        'count_matrix_over_aggregate.mtx',
        'count_matrix_over_aggregate.rownames',
        'count_matrix_over_aggregate.colnames',
        'qc_metrics/dup_level.txt',
        'qc_metrics/mapping_rate.txt',
        'qc_metrics/mt_content.txt',
        'qc_metrics/sequencing_depth.txt',
        'qc_metrics/uniq_nuc_frags.txt',
        'qc_metrics/frip.txt',
        'qc_metrics/frac_open.txt',
        'qc_metrics/library_size.txt'

rule fastp:
    input:
        r1='{rep}/fastq/{cell}_r1.fq.gz',
        r2='{rep}/fastq/{cell}_r2.fq.gz'
    output:
        r1='{rep}/fastp_fq/{cell}_r1_trimmed.fq.gz',
        r2='{rep}/fastp_fq/{cell}_r2_trimmed.fq.gz',
        h='{rep}/fastp_report/{cell}_fastp.html',
        j='{rep}/fastp_report/{cell}_fastp.json'
    log:
        out='logs/fastp/{rep}/{cell}.stdout',
        err='logs/fastp/{rep}/{cell}.stderr'
    threads: 12
    shell:
        ''' fastp -l 25 -w {threads} --detect_adapter_for_pe \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.h} -j {output.j} \
            1> {log.out} 2> {log.err}
        '''

rule hisat2:
    input:
        r1='{rep}/fastp_fq/{cell}_r1_trimmed.fq.gz',
        r2='{rep}/fastp_fq/{cell}_r2_trimmed.fq.gz'
    output:
        bam='{rep}/hisat2_mapped/{cell}_f2q30.bam',
        stats='{rep}/hisat2_log/{cell}_aln_sum.txt'
    params:
        idx=config['genome'],
        maxi=config['hisat2_X']
    threads: 12
    shell:
        ''' hisat2 \
            -X {params.maxi} \
            -p {threads} \
            --no-temp-splicesite \
            --no-spliced-alignment \
            -x {params.idx} \
            -1 {input.r1} \
            -2 {input.r2} \
            -3 1 \
            --summary-file {output.stats} | \
            samtools view -ShuF 4 -f 2 -q 30 - | \
            samtools sort - -T {wildcards.cell}_tmp -o {output.bam}
        '''

rule spicard:
    input:
        bam='{rep}/hisat2_mapped/{cell}_f2q30.bam',
        pmd=config['picard_jar']
    output:
        bam='{rep}/picard_bam/{cell}_f2q30_pmd.bam',
        met='{rep}/picard_log/{cell}_f2q30_pmd.out',
    log:
        out='logs/spicard/{rep}/{cell}.stdout',
        err='logs/spicard/{rep}/{cell}.stderr'
    shell:
        ''' java -jar -Xmx4g {input.pmd} \
            MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            METRICS_FILE={output.met} \
            1> {log.out} 2> {log.err}
        '''

rule index:
    input:
        '{rep}/picard_bam/{cell}_f2q30_pmd.bam'
    output:
        '{rep}/picard_bam/{cell}_f2q30_pmd.bam.bai'
    shell:
        ''' samtools index {input}
        '''

rule isize:
    input:
        '{rep}/picard_bam/{cell}_f2q30_pmd.bam'
    output:
        '{rep}/isize_hist/{cell}_isize.hist'
    shell:
        """ samtools view {input} | \
            sed '/chrM/d' | \
            awk '{{ if($9>0){{print $9}} else {{print -1*$9}} }}' | \
            sort | uniq -c | \
            sort -b -k2,2n | \
            sed -e 's/^[ \t]*//' > {output}
        """

rule list_bam:
    input:
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam', zip,
               rep=samples["rep"],
               cell=samples["cell"])
    output:
        expand('{rep}/bam_file_list.txt', rep=replicate)
    shell:
        ''' scripts/list_bam.sh
        '''

# there is a limit of maximum number of files you can
# open at once (ulimit -n, which is often 1024)
# if you have many cells, better merge them gradually

rule merge_rep:
    input:
        '{rep}/bam_file_list.txt'
    output:
        '{rep}/f2q30_merged.bam'
    shell:
        ''' samtools merge -b {input} {output}
        '''

rule merge_all:
    input:
        expand('{rep}/f2q30_merged.bam', rep=replicate)
    output:
        'aggregate/f2q30_merged.bam'
    shell:
        ''' samtools merge {output} {input}
        '''

rule mpicard:
    input:
        bam='aggregate/f2q30_merged.bam',
        pmd=config['picard_jar']
    output:
        bam='aggregate/f2q30_merged_pmd.bam',
        met='aggregate/f2q30_merged_pmd.out'
    log:
        out='logs/mpicard/mpicard.stdout',
        err='logs/mpicard/mpicard.stderr'
    shell:
        ''' java -jar -Xmx8g {input.pmd} \
            MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            METRICS_FILE={output.met} \
            1> {log.out} 2> {log.err}
        '''

rule isize_aggregate:
    input:
        'aggregate/f2q30_merged_pmd.bam'
    output:
        'aggregate/f2q30_merged_pmd_isize.hist'
    shell:
        """ samtools view {input} | \
            sed '/chrM/d' | \
            awk '{{ if($9>0){{print $9}} else {{print -1*$9}} }}' | \
            sort | uniq -c | \
            sort -b -k2,2n | \
            sed -e 's/^[ \t]*//' > {output}
        """

rule bam2bed:
    input:
        bam='aggregate/f2q30_merged_pmd.bam',
        bl=config['blacklist']
    output:
        'aggregate/f2q30_merged_pmd.bed'
    shell:
        ''' bamToBed -i {input.bam} | \
            intersectBed -a - -b {input.bl} -v > {output}
        '''

rule macs2:
    input:
        'aggregate/f2q30_merged_pmd.bed'
    output:
        'aggregate/aggregated_scATAC_peaks.narrowPeak',
        'aggregate/aggregated_scATAC_peaks.xls',
        'aggregate/aggregated_scATAC_treat_pileup.bdg',
        'aggregate/aggregated_scATAC_control_lambda.bdg'
    params:
        gs=config['gsize'],
        broad=config['bpk'],
        fmt=config['macs2_format'],
        shift=config['macs2_shift']
    log:
        out='logs/macs2/macs2.stdout',
        err='logs/macs2/macs2.stderr'
    shell:
        ''' macs2 callpeak -t {input} \
            -g {params.gs} \
            {params.broad} \
            -f {params.fmt} \
            -q 0.01 \
            --nomodel \
            {params.shift} \
            --keep-dup all \
            -B --SPMR \
            --outdir aggregate \
            -n aggregated_scATAC \
            1> {log.out} 2> {log.err}
        '''

rule bigwig:
    input:
        'aggregate/aggregated_scATAC_treat_pileup.bdg',
        config['chromsize']
    output:
        'aggregate/aggregated_scATAC_treat_pileup.bw'
    shell:
        ''' bdg2bw {input}
        '''

rule count:
    input:
        peak='aggregate/aggregated_scATAC_peaks.narrowPeak',
        bam='{rep}/picard_bam/{cell}_f2q30_pmd.bam'
    output:
        '{rep}/count/{cell}.count'
    shell:
        ''' coverageBed \
            -a {input.peak} \
            -b {input.bam} | \
            cut -f 4,11 > {output}
        '''

rule countMatrix:
    input:
        expand('{rep}/count/{cell}.count', zip,
               rep=samples["rep"], cell=samples["cell"])
    output:
        'count_matrix_over_aggregate.mtx',
        'count_matrix_over_aggregate.rownames',
        'count_matrix_over_aggregate.colnames'
    script:
        'scripts/genernate_count_matrix.py'

rule basicQc:
    input:
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam', zip,
               rep=samples["rep"], cell=samples["cell"]),
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam.bai', zip,
               rep=samples["rep"], cell=samples["cell"]),
        expand('{rep}/hisat2_log/{cell}_aln_sum.txt', zip,
               rep=samples["rep"], cell=samples["cell"])
    output:
        'qc_metrics/dup_level.txt',
        'qc_metrics/mapping_rate.txt',
        'qc_metrics/mt_content.txt',
        'qc_metrics/sequencing_depth.txt',
        'qc_metrics/uniq_nuc_frags.txt',
        'qc_metrics/library_size.txt'
    shell:
        ''' scripts/get_dup_level.sh
            scripts/get_depth_mr.sh
            scripts/get_ufrags_mt.sh
            scripts/get_lib_size.sh
        '''

rule frip:
    input:
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam', zip,
               rep=samples["rep"], cell=samples["cell"]),
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam.bai', zip,
               rep=samples["rep"], cell=samples["cell"]),
        'aggregate/aggregated_scATAC_peaks.narrowPeak'
    output:
        'qc_metrics/frip.txt'
    shell:
        ''' scripts/get_frip.sh
        '''

rule fracOpen:
    input:
        expand('{rep}/picard_bam/{cell}_f2q30_pmd.bam', zip,
               rep=samples["rep"], cell=samples["cell"]),
        'aggregate/aggregated_scATAC_peaks.narrowPeak'
    output:
        'qc_metrics/frac_open.txt'
    shell:
        ''' scripts/get_frac_open.sh
        '''
