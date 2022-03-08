#!/bin/bash
​
source ~/.bashrc
conda activate merida
​
snakemake \
    --latency-wait 100 \
    --jobs 1000 \
    --cluster "sbatch -J {params.jobname} -t {params.runtime} -c {params.cpu} --mem={params.mem} --output {params.slurm_output} --mail-type {cluster.mail_type} --mail-user {cluster.email}" \
    --cluster-config scripts/slurm/cluster_config.yaml
