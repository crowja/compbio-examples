#! /bin/bash

. "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bioconductor-env

./ex1.R

