import pandas as pd
from glob import glob

files = glob('qc_metrics/*.txt')

dfs = []

for f in files:
    dfs.append(pd.read_csv(f, sep='\t', index_col=0))

out = pd.concat(dfs, axis=1)
out.to_csv('outs/sample_info.csv')
