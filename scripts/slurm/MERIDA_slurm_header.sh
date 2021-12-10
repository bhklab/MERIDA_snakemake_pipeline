#!/bin/bash
#SBATCH --account=def-bhaibeka-ab
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --job-name=full_MERIDA
#SBATCH --output=full_MERIDA_%j.out
/home/emilyso/scratch/MERIDA/build/MERIDA_ILP /home/emilyso/scratch/MERIDA/Example_Config1.txt no
