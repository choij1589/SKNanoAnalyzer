#!/bin/bash
export PATH="/opt/conda/bin:${PATH}"
export MAMBA_ROOT_PREFIX="/opt/conda"
eval "$(micromamba shell hook -s bash)"
micromamba activate Nano

cd [WORKDIR]
hadd -f -j 8 [TARGET] output/hists_*.root
exit $?
