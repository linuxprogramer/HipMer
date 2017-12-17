#!/bin/bash
# Set SGE options:
#$ -cwd
#$ -j y
#$ -pe pe_16 16
#$ -l ram.c=7.5G

export USE_SHM=1
export UPC_SHARED_HEAP_MB=1500

BASE=$(dirname $(which $0))
source ${BASE}/test_hipmer_chr14-genepool.sh
