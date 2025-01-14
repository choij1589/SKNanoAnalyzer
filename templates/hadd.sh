#!/bin/bash
source /opt/conda/bin/activate
conda activate torch

cd [WORKDIR]
hadd -f -j 8 [TARGET] output/hists_*.root
exit $?
