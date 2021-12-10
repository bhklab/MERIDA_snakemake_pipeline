#! /bin/bash 
​
source ~/.bashrc
conda activate merida
​
snakemake \
    --latency-wait 100 \
    --jobs 30 \
    --cluster "sbatch -J {params.jobname} -t {params.runtime} --mail-type {cluster.mail_type} --mail-user {cluster.email}" \
    --cluster-config scripts/slurm/cluster_config.yaml