#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

export USE_SHM=1
export UPC_SHARED_HEAP_MB=1000

source .edison_deploy/test_hipmer_chr14-edison.sh
