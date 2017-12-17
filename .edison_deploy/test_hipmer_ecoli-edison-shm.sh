#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

export USE_SHM=1

source .edison_deploy/test_hipmer_ecoli-edison.sh


