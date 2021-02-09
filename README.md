# A Snakemake pipeline to process plate-based scATAC-seq data
This repository contains codes for processing scATAC-seq data produced by the [plate-based scATAC-seq method](https://www.nature.com/articles/s41467-018-07771-0).

## What are the differences comparing to the original Nat. Comms. publication?

In this pipeline, we have:

1. Added `config.json` file to make the processing more flexible and easier to modify.
2. Used the `BED` file for the MACS2 peak calling. See [this tweet](https://twitter.com/XiChenUoM/status/1336658454866325506) for the reason.
3. Produced `10x Genomics` like output files in the final `outs` directory so that they can be easily put into downwstream analysis package like [Signac](https://satijalab.org/signac/). See below for more details.

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

Get `calc`, `addCols`, `bedClip`, `bedGraphToBigWig` and `fetchChromSizes` from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/), and make sure they are executable and in your `$PATH`.

Get [bdg2bw](https://gist.github.com/taoliu/2469050) to convert the macs2 generated begraph to bigwig for visualisation.

### 2. Change the content insde `config.json`

There are a few things you need to change in the `config.json` file according to your computing environment:

__genome__: This is the prefix of `hisat2` index for the genome you are working on. It is basically passed to the `-x` flag of `hisat2` during alignment.

__hisat2_X__: This is the `-X` flag of `hisat2`, which suggests the maximum frament in bp allowed during the alignment. People normally use `2000` for ATAC-seq.

__picard_jar__: The location of the `picard.jar` file.

__blacklist__: The ENCODE blacklist region to exclude for analysis. Check [this publication](https://www.nature.com/articles/s41598-019-45839-z) for more details. There are a few pre-compiled blacklists for different genome builds that can be found [here](https://github.com/Boyle-Lab/Blacklist).

__gsize__: This is the genome size `-g` option for `macs2` during the peak calling.

__bpk__: This contains the `macs2` flags for broadPeak calling. In most cases, ATAC-seq signals are sharp, we normally leave this empty here.

__chromsize__: The location pointing to the tab-delimited file that contains the length of each chromosome. Use the UCSC `fetchChromSizes` program to get this file. For example, to get the file for hg38, simply run `fetchChromSizes hg38 > hg38.chrom.sizes`.

__macs2_format__: This is the file format duing `macs2` peak calling. We use `BED` in this pipeline.

__macs2_shift__: The flags used for calling narrowPeak. Use `--nomodel --shift -100 --extsize 200` to centre the reads on the Tn5 cutting sites.

### 3. Organise your files

The starting point of the pipeline is the `fastq` files. Put your `fastq` files inside each plate directory. Also put the `Snakefile`, `config.json` and the `scripts` folder from this repository to your experiment directory. The structure will be like this:

```
Experiment
│
├── config.json
│
├── cluster.json
│
├── submit_snake.sh
│
├── plate1
│   │
│   ├── fastq
│   │   │
│   │   ├── scATAC_p1_cell_001_r1.fq.gz
│   │   ├── scATAC_p1_cell_001_r2.fq.gz
│   │   ├── scATAC_p1_cell_002_r1.fq.gz
│   │   ├── scATAC_p1_cell_002_r2.fq.gz
│   │   ├── scATAC_p1_cell_003_r1.fq.gz
│   │   ├── scATAC_p1_cell_003_r2.fq.gz
.   .   .
.   .   .
├── plate2
│   │
│   ├── fastq
│   │   │
│   │   ├── scATAC_p2_cell_001_r1.fq.gz
│   │   ├── scATAC_p2_cell_001_r2.fq.gz
│   │   ├── scATAC_p2_cell_002_r1.fq.gz
│   │   ├── scATAC_p2_cell_002_r2.fq.gz
│   │   ├── scATAC_p2_cell_003_r1.fq.gz
│   │   ├── scATAC_p2_cell_003_r2.fq.gz
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

### Understanding the output files

There will be quite a few useful intermediate files generated during the process of the pipeline. They are organised into each directory, and the name of the file is self-explanatory. The most important output files are in the `outs` directory. If the the pipeline runs successfully, you should expect an `outs` directory under the `Experiment` directory. Inside the `outs` directory, there will be six files:

__aggregate_fragments.tsv.gz__: this is a tab-delimited file that contains the ATAC fragments of all cells after deduplication, with the following specification:

| column     | meaning                                     |
|------------|---------------------------------------------|
| 1st column | chromosome of the fragment                  |
| 2nd column | 0-based start coordinate of the fragment    |
| 3rd column | 1-based end corrdinate of the fragment      |
| 4th column | the cell name from where the fragment comes |
| 5th column | just '1'                                    |

__aggregate_fragments.tsv.gz.tbi__: the index of the fragment file, created by the `indexFrag` rule from the `Snakefile`.

__count_matrix_over_aggregate.mtx__: the peak-by-cell count matrix in `matrix market format`. This is basically `sparse.csc_matrix` if you use `python`; or this can be treated as `dgCMatrix` if you use `R`.

__count_matrix_over_aggregate.cols__: the name of each cell in plain text format.

__count_matrix_over_aggregate.rows__: the peak location in a 3-column `bed` format.

__sample_info.csv__: a csv file containing the basic quality metrics of each cell. The meaning of each column is described as follows:

| column           | value                                                                                      | typical range for a successful cell |
|------------------|--------------------------------------------------------------------------------------------|-------------------------------------|
| cell             | the name of the cell                                                                       | N/A                                 |
| frac_open        | percentage (%) of all peaks detected (at least one read) in the cell                       | 1 - 20                              |
| mapping_rate     | overall alignment rate (%) from hisat2                                                     | 70 - 99                             |
| mt_content       | percentage (%) of reads mapped to the reference genome                                     | 0.1 - 90                            |
| uniq_nuc_frags   | number of read mapped to the nuclear genome after deduplication                            | 10,000 - 100,000                    |
| dup_level        | duplication level estimated by the picard tool, indicating the fraction of duplicate reads | 0.4 - 0.9                           |
| frip             | percentage (%) of reads that come from the peak region                                     | 20 - 80                             |
| sequencing_depth | total number of reads sequenced per cell                                                   | 10,000 - 1,000,000                  |
| library_size     | library complexity (number of unique fragments) estimated by the picard tool               | 10,000 - 1,000,000                  |

### 5. Load the output files into [Signac](https://satijalab.org/signac/)

If you use `python`, use `mmread` from `scipy` to load the `mtx` file and conduct analysis using differen packages from `scikit-learn`. You can also try [EpiScanpy](https://episcanpy.readthedocs.io/en/latest/index.html). If you use `R`, you have many choices for the analysis. To load data into `Signac`, use the following lines of code:

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(dplyr)
library(readr)

# read the content from the 'outs' directory
setwd("/your/working/directory")
mex_dir_path <- "/path/to/mtx"

mtx_path <- paste(mex_dir_path, "count_matrix_over_aggregate.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "count_matrix_over_aggregate.rows", sep = '/')
barcode_path <- paste(mex_dir_path, "count_matrix_over_aggregate.cols", sep = '/')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
metadata <- read.csv(
  file = "/path/to/outs/sample_info.csv",
  header = TRUE,
  row.names = 1
)

# create a Signac chromatin assay and a Seurat object
mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)
chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c("_", "_"),
  genome = 'hg38',
  min.cells = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  project = 'scATAC-seq_is_cool',
  meta.data = metadata
)
```

# Contact
Xi Chen  
chenx9@sustech.edu.cn
