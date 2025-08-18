#!/bin/bash
RUN=$1

SKNano.py -a ClosDiLepTrigs -i DYJets,TTLL_powheg -n 30 -r $RUN --userflags RunDiMu
SKNano.py -a ClosDiLepTrigs -i DYJets,TTLL_powheg -n 30 -r $RUN --userflags RunEMu
if [[ $RUN == "Run2" ]]; then
    SKNano.py -a ClosDiLepTrigs -i WZTo3LNu_amcatnlo,TTZToLLNuNu -n 30 -r $RUN --userflags Run1E2Mu
    SKNano.py -a ClosDiLepTrigs -i WZTo3LNu_amcatnlo,TTZToLLNuNu -n 30 -r $RUN --userflags Run3Mu
elif [[ $RUN == "Run3" ]]; then
    SKNano.py -a ClosDiLepTrigs -i WZTo3LNu_powheg,TTZ_M50,TTZ_M4to50 -n 30 -r $RUN --userflags Run1E2Mu
    SKNano.py -a ClosDiLepTrigs -i WZTo3LNu_powheg,TTZ_M50,TTZ_M4to50 -n 30 -r $RUN --userflags Run3Mu
else
    echo "Unknown run: $RUN"
    exit 1
fi
