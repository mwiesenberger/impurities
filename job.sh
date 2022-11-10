#!/bin/bash

#SBATCH -N 1 -n 1 --ntasks-per-node=1
#SBATCH --gres=gpu:1
####SBATCH --gpus-per-node=1
#SBATCH --partition=m100_fua_prod
#SBATCH --account=fuac6_tsvv3
#SBATCH --time=24:00:00 # 24 hours is maximum

echo "Marconi 100 cluster with Tesla V100 GPUs"

hostname
date
module list
echo "$@"

# $@ forwards all arguments
./impurities_hpc "$@"

date
