#! /bin/bash

. "$(conda info --base)/etc/profile.d/conda.sh"
##conda activate bioconductor-env
conda activate deseq2-env

./ex1.R

