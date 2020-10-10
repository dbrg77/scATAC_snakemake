snakemake --timestamp -j 500 --cluster-config cluster.json --cluster "bsub -M {cluster.memory} -n {cluster.nCPUs} -R {cluster.resources} -o {cluster.output}"
