# Plate-based scATAC-seq vs 10x Genomics Single Cell ATAC
This directory contains information to reproduce the comparison results in the paper.

### 1. Downsample the plate-based scATAC-seq data

[seqtk](https://github.com/lh3/seqtk) was used to downsample the plate-based scATAC-seq data from human PBMC to around 50k reads per cell, which is around 40% of the full data. The following commands were used:

```
seqtk sample -s 42 {cell}_r1.fq.gz 0.4 | gzip > {cell}_ds_r1.fq.gz
seqtk sample -s 42 {cell}_r2.fq.gz 0.4 | gzip > {cell}_ds_r2.fq.gz
```

Then, the downsampled fastq files can be used as the starting point for the snakemake pipeline.

### 2. Get original fastq files from the 10x Genomics website

There are quite a few scATAC-seq datasets from the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/datasets), we chose the two 500-cell datasets using the v1 and the NextGEM chemistries, as these two are the best among all the PBMC datasets.

To get [the data from the v1 chemistry](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_pbmc_500_v1), download the fastq and the metadata:

```
wget -c https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_500_v1/atac_pbmc_500_v1_fastqs.tar

wget -c https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_500_v1/atac_pbmc_500_v1_singlecell.csv
```

Extract the fastq tar, you will get the following four files:

__atac_pbmc_500_v1_S1_L001_I1_001.fastq.gz__: this is the sample index file, which can be discarded.  
__atac_pbmc_500_v1_S1_L001_R1_001.fastq.gz__: this is the Read 1 file of ATAC-seq.  
__atac_pbmc_500_v1_S1_L001_R2_001.fastq.gz__: this is the 10x cell barcodes read file.  
__atac_pbmc_500_v1_S1_L001_R3_001.fastq.gz__: this is the Read 2 file of ATAC-seq.

Creat an `index.txt` file which contains the "true cell barcodes" using the information from the metadata csv file:

```
echo -e "#index1\tName" > index.txt
awk -F, '$10==1' atac_pbmc_500_v1_singlecell.csv | cut -f 1 -d, | cut -f 1 -d '-' | awk 'BEGIN{OFS="\t"}{print $1, "cell_" NR}' >> index.txt
```

Demultiplex the fastq files into individual true cells using [deML](https://github.com/grenaud/deml):

```
deML -i index.txt -f atac_pbmc_500_v1_S1_L001_R1_001.fastq.gz -r atac_pbmc_500_v1_S1_L001_R3_001.fastq.gz -if1 atac_pbmc_500_v1_S1_L001_R2_001.fastq.gz --mm 1 -o pbmc_v1
```

After the above step, remove all 'failed' and 'unknown' files. You will get two fastq files per cell: `{cell}_r1.fq.gz` and `{cell}_r2.fq.gz`. These can be used as the starting point for the snakemake pipeline.

To get [the data from the NextGEM chemistry](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_pbmc_500_nextgem), download the fastq and the metadata:

```
wget -c https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fastqs.tar

wget -c https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_singlecell.csv
```

Then you use the exact the same procedures used for the v1 data to demultiplex fastq files into each individual true cell barcodes and run the snakemake pipeline.

### 3. Collect the quality metrics

In each experiment, there is a `sample_info.csv` in their own `outs` directory. The metrics in those files are used to make the plot in the paper.

# Contact
Xi Chen  
chenx9@sustech.edu.cn
