import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import array
import pysam
import sys

bam = pysam.AlignmentFile(sys.argv[1], "rb")
pdf = ''.join([sys.argv[2], "_isize_hist.pdf"])
txt = ''.join([sys.argv[2], "_isize_density_xy_values.txt"])

frags = array.array('i', [])

for read in bam:
    
    if (not read.is_unmapped and \
        read.reference_name != "chrM" and \
        read.is_proper_pair and \
        read.is_read1 > 0):
        
        frags.append(abs(read.template_length))

bam.close()

# plot histogram
plt.subplots(figsize=(8,8))
y, x, p = plt.hist(frags, bins=range(0,2001), range=(0,2000),
                   normed=1, color='r', histtype="step", alpha=.8)
plt.xticks(range(0,2500,500), map(str, range(0,2500,500)), fontsize=16)
plt.xlabel("Length (bp)", fontsize=18)
plt.yticks(fontsize=16)
plt.ylabel("Density", fontsize=18)
plt.tight_layout()
plt.savefig(pdf, transparent=True)

# output x,y values for more flexible plot and qc by other programs
with open(txt, 'w') as fh:
    fh.write("length\tdensity\n")
    for i,j in enumerate(y):
        fh.write("%s\t%s\n" % (x[i], j))
