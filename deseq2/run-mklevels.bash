#! /bin/bash

. "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scicomp

./mklevels.py

