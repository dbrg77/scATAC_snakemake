###################################################################################
# Xi Chen (xi.chen.xchen at gmail dot com)
#
# usage:
#     python Fragment_length_density_plot.py <bam_files> <labels> <out_prefix>
# 
# output:
#     out_prefix_histogram.pdf
#     and
#     out_prefix_log_scale.pdf
#
# example:
#     python Fragment_length_density_plot.py ESC.bam NPC.bam ESC NPC density_plot
# 
# it will generate two files:
# density_plot_histogram.png and density_plot_log_scale.png
# 
# The fragment length info. from the two bam input files
#   will be plotted in the same figure.
#
#
##################################################################################

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

d = len(sys.argv)/2

files = sys.argv[1:d]
labels = sys.argv[d:-1]
prefix = sys.argv[-1]
colour = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

bins = np.linspace(0, 1200, 1201)
data = np.empty((len(files), 2, 1200))

frags = []

for i,j in enumerate(files):
    frags.append([])
    bam = pysam.Samfile(j, 'rb')
    for read in bam:
        if bam.getrname(read.tid)!="chrM" and read.tlen>50 and read.tlen<1500:
            frags[i].append(read.tlen)
    bam.close()
    hist, b = np.histogram(frags[i], bins, density=True)
    data[i][0] = (b[1:]+b[:-1])/2
    data[i][1] = hist

for i,j in enumerate(data):
    xs, ys = j
    plt.plot(xs, ys, colour[i]+'-', label=labels[i], alpha=0.7)

plt.xticks(range(0,1600,200))
plt.xlim(0,1600)
plt.xlabel('Fragment length (bp)')
plt.ylabel('Density')
plt.legend()
plt.savefig(prefix+"_histogram.pdf")

plt.clf()
ymax = 0

for i,j in enumerate(data):
    xs, ys = j
    if max(ys) > ymax:
        ymax = max(ys)
    plt.plot(xs, ys, colour[i]+'-', label=labels[i], alpha=0.7)

plt.xticks(range(0,1600,200))
plt.xlim(0,1600)
plt.yscale('log')
plt.ylim(ymax/10000., ymax+0.05)
plt.xlabel('Fragment length (bp)')
plt.ylabel('Density')
plt.legend()
plt.savefig(prefix+"_log_scale.pdf")
