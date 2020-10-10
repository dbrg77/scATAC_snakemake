snakemake -j 12 --cluster-config cluster.json --cluster "bsub -q {cluster.queues} -M {cluster.memory} -n {cluster.nCPUs} -R {cluster.resources} -o {cluster.output}"
