# ATAC

Scripts in this repos are used to perform ATAC-seq analyses routines.
```sh
frag_distr.py
```

# Requirements
  - matplotlib
  - pysam

# Usage
```sh
python frag_distr.py input_alignment.bam out_prfix
```

# Output:

* out_prefix_isize_hist.pdf: density plot of insert size.
* out_prefix_isize_density_xy_values.txt: a two-column tab-delimited file with length:density pair.
