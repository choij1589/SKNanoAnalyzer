#!/bin/bash
CHANNEL=$1

## Update DATASTREAM
DATASTREAM="Skim_TriLep_MuonEG"
MC_LIST="SampleLists/Run2NanoV9"

SKNano.py -a PromptAnalyzer -i ${DATASTREAM} -n 10 -e 2018 --userflags ${CHANNEL},RunNoHEMVeto,NoTreeMode --python --memory 4000
SKNano.py -a MatrixAnalyzer -i ${DATASTREAM} -n 10 -e 2018 --userflags ${CHANNEL},RunNoHEMVeto,NoTreeMode --python --memory 4000
SKNano.py -a PromptAnalyzer -i ${MC_LIST}/TriLepton.txt -n 40 -e 2018 --userflags ${CHANNEL},RunSyst,RunNoHEMVeto,NoTreeMode --python --memory 4000
