# output csc_matrix format
## rows are peaks, columns are cells

import subprocess
from glob import iglob

scs = iglob('*/count/*.count')

header = '%%MatrixMarket matrix coordinate integer general\n'
spacer = '%\n'

mtx = open('outs/count_matrix_over_aggregate.mtx', 'w')
mtx.write(header)
mtx.write(spacer)

rownames = open('outs/count_matrix_over_aggregate.rows', 'w')
colnames = open('outs/count_matrix_over_aggregate.cols', 'w')

with open('aggregate/aggregated_scATAC_peaks_formatted.bed') as peak:
    for line in peak:
        pos = line.split('\t')[:3]
        rownames.write('\t'.join(pos))
        rownames.write('\n')

total_non_zeroes = 0

for c, j in enumerate(scs):
    pl = j.split('/')[0]
    fn = j.split('/')[-1]
    colnames.write(pl + '_' + fn[:-6] + '\n')
    with open(j) as fh:
        for p, v in enumerate(fh):
            pid, count = v.strip().split('\t')
            coln = c + 1
            if count != '0':
                total_non_zeroes += 1
                rown = p + 1
                mtx.write('%s %s %s\n' %(rown, coln, count))

total_c = c + 1
total_p = p + 1

mtx.close()
rownames.close()
colnames.close()

# add the information of numbers of rows, columns and non-zeroes
# to the third (the "3i" command of sed) line of the mtx file
subprocess.run("sed -i '3i%s %s %s' outs/count_matrix_over_aggregate.mtx" %
               (total_p, total_c, total_non_zeroes),
               shell=True)
