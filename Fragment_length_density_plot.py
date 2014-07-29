#############################################################################
#
#
# usage:
#     python Fragment_length_density_plot <bam_files> <labels> <out_prefix>
# output:
#     out_prefix_histogram.png
#     and
#     out_prefix_log_scale.png
#
#
#############################################################################

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
        if bam.getrname(read.tid)!="chrM" and read.tlen>0 and read.tlen<1200:
            frags[i].append(read.tlen)
    bam.close()
    hist, b = np.histogram(frags[i], bins, density=True)
    data[i][0] = (b[1:]+b[:-1])/2
    data[i][1] = hist

for i,j in enumerate(data):
    xs, ys = j
    plt.plot(xs, ys, colour[i]+'-', label=labels[i], alpha=0.7)

plt.xticks(range(0,1200,200))
plt.xlim(0,1200)
plt.xlabel('Fragment length (bp)')
plt.ylabel('Density')
plt.legend()
plt.savefig(prefix+"_histogram.png")

plt.clf()
ymax = 0

for i,j in enumerate(data):
    xs, ys = j
    if max(ys) > ymax:
        ymax = max(ys)
    plt.plot(xs, ys, colour[i]+'-', label=labels[i], alpha=0.7)

plt.xticks(range(0,1200,200))
plt.xlim(0,1200)
plt.yscale('log')
plt.ylim(ymax/10000., ymax+0.05)
plt.xlabel('Fragment length (bp)')
plt.ylabel('Density')
plt.legend()
plt.savefig(prefix+"_log_scale.png")
