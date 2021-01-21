# A Snakemake pipeline to process plate-based scATAC-seq data
This repository contains codes for processing scATAC-seq data produced by the [plate-based scATAC-seq method](https://www.nature.com/articles/s41467-018-07771-0).

## How to use?

### 1. Get all the required softwares/packages

Install the following packages either using conda or pip or directly download the pre-compiled binary from the website:

```
python3
snakemake 5.3.0
numpy v1.18.5
scipy v1.5.0
pandas v1.0.5
hisat2 v2.1.0
samtools v1.9
bedtools v2.27.1
fastp v0.20.1
macs2 v2.2.7.1
tabix 0.2.5
```

Get the picard tool `picard.jar` from https://github.com/broadinstitute/picard/releases

Get `calc`, `addCols`, `bedClip`, `bedGraphToBigWig` and `fetchChromSizes` from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/), and make sure they are executable and in your `PATH`.

Get [bdg2bw](https://gist.github.com/taoliu/2469050) to convert the macs2 generated begraph to bigwig for visualisation.

### 2. Change the content insde `config.json`

There are a few things you need to change in the `config.json` according to your computing environment:

__genome__: This is the prefix of `hisat2` index for the genome you are working on. It is basically passed to the `-x` flag of `hisat2` during alignment.

__hisat2_X__: This is the `-X` of `hisat2`, which suggest the maximum frament in bp allowed during the alignment. People normally use `2000` for ATAC-seq.

__picard_jar__: The location pointing to the `picard.jar` file is.

__blacklist__: The ENCODE blacklist region to exclude for analysis. Check [this publication](https://www.nature.com/articles/s41598-019-45839-z) for more details. There are a few pre-compiled blacklists for different genome builds that can be found [here](https://github.com/Boyle-Lab/Blacklist).

__gsize__: This is the `-g` genome size option for `macs2` during peak calling.

__bpk__: This contains the `macs2` flags for broadPeak calling. In most cases, ATAC-seq signals are sharp, we normally leave this empty here.

__chromsize__: The location pointing to the tab-delimited file that contains the length of each chromosome. Use the UCSC `fetchChromSizes` to get this file. For example, to get the file for hg38, simply run `fetchChromSizes hg38 > hg38.chrom.sizes`.

__macs2_format__: This is the file format duing `macs2` peak calling. We use `BED` in this pipeline.

__macs2_shift__: The flags used for calling narrowPeak. Use `--nomodel --shift -100 --extsize 200` to centre the reads on the Tn5 cutting sites.

### 3. Organise your files

The starting point of the pipeline is the `fastq` files. Put your `fastq` files inside each plate directory. Also put the `Snakefile`, `config.json` and the `scripts` folder from this repositor to your experiment directory. The structure will be like this:

```
Experiment
│
├── config.json
│
├── plate1
│   │
│   ├── fastq
│   │   │
│   │   ├── scATAC_cell_001_r1.fq.gz
│   │   ├── scATAC_cell_001_r2.fq.gz
│   │   ├── scATAC_cell_002_r1.fq.gz
│   │   ├── scATAC_cell_002_r2.fq.gz
│   │   ├── scATAC_cell_003_r1.fq.gz
│   │   ├── scATAC_cell_003_r2.fq.gz
.   .   .
.   .   .
├── plate2
│   │
│   ├── fastq
│   │   │
│   │   ├── scATAC_cell_001_r1.fq.gz
│   │   ├── scATAC_cell_001_r2.fq.gz
│   │   ├── scATAC_cell_002_r1.fq.gz
│   │   ├── scATAC_cell_002_r2.fq.gz
│   │   ├── scATAC_cell_003_r1.fq.gz
│   │   ├── scATAC_cell_003_r2.fq.gz
.   .   .
.   .   .
├── scripts
│   │
│   ├── collect_metadata.py
│   ├── generate_count_csc_mtx.py
│   ├── generate_fragments_file.sh
│   ├── get_depth_mr.sh
│   ├── get_dup_level.sh
│   ├── get_frac_open.sh
│   ├── get_frip.sh
│   ├── get_lib_size.sh
│   ├── get_ufrags_mt.sh
│   └── list_bam.sh
│
└── Snakefile

```

### 4. Run the processing pipeline

To use all available cores to run the pipeline, simply type `snakemake --cores` under the `Experiment` directory.

You can also run the pipeline using `bsub`, using the command and setting provided in the `submit_snake.sh` and `cluster.json` files. These two files can be ignored if you are not using `bsub`.

## What are the differences?

There are not many differences comparing to the [original analysis method](https://github.com/dbrg77/plate_scATAC-seq) used in the original publication. In the original repository, the data processing is hard coded for a specific species. Here we have added `config.json` file to make the processing more flexible and easier to modify.

# Contact
Xi Chen  
chenx9@sustech.edu.cn
