ATAC
====

Scripts in this repos are used to perform ATAC-seq analyses routines.

Fragment_length_density_plot.py
-------------------------------

Requirements:
~~~~~~~~~~~~~~

matplotlib
numpy
pysam
~~~~~~~~~~~~~

Usage:

python Fragment_length_density_plot.py <input_bam_files> <sample_labels> <output_file_prefix>

The script takes multiple paired-end bam files as input, extract the isize information from 
properly paired reads, and plot the density of the isize. It can plot a maximum of 7 samples
in one figure.

Example:


Plot one sample:

python Fragment_length_density_plot.py input1.bam sample1 density

The command generates two figures: density_hitogram.png (regular scale) and density_log_scale.png (log scale)

Plot three samples:

python Fragment_length_density_plot.py input1.bam input2.bam input3.bam sample1 sample2 sample3 out

The command generates two figures: out_hitogram.png and out_log_scale.png. Each figures contains three samples (1 - 3).
