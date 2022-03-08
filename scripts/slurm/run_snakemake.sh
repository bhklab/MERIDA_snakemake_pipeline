#! /bin/bash
​
source ~/.bashrc
conda activate merida
​
snakemake -n -s Snakefile \
    --latency-wait 100 \
    -j 30
    --cluster "sbatch -J {params.jobname} -t {params.runtime} --output {params.slurm_output} --mail-type {cluster.mail_type} --mail-user {cluster.email}"
    --cluster-config scripts/slurm/cluster_config.yaml