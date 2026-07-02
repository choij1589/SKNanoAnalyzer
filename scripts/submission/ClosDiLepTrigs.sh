#!/bin/bash
RUN=$1

SKNano.py -a ClosDiLepTrigs -i DYJets,TTLL_powheg -n 30 -r $RUN --userflags RunDiMu
SKNano.py -a ClosDiLepTrigs -i DYJets,TTLL_powheg -n 30 -r $RUN --userflags RunEMu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc70_MA15 -n 10 -r $RUN --userflags Run1E2Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc100_MA60 -n 10 -r $RUN --userflags Run1E2Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc130_MA90 -n 10 -r $RUN --userflags Run1E2Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc160_MA155 -n 10 -r $RUN --userflags Run1E2Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc70_MA15 -n 10 -r $RUN --userflags Run3Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc100_MA60 -n 10 -r $RUN --userflags Run3Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc130_MA90 -n 10 -r $RUN --userflags Run3Mu
SKNano.py -a ClosDiLepTrigs -i TTToHcToWAToMuMu-MHc160_MA155 -n 10 -r $RUN --userflags Run3Mu


if [[ $RUN == "Run2" ]]; then
    SKNano.py -a ClosDiLepTrigs -i Skim_TriLep_WZTo3LNu_amcatnlo,Skim_TriLep_TTZToLLNuNu -n 10 -r $RUN --userflags Run1E2Mu
    SKNano.py -a ClosDiLepTrigs -i Skim_TriLep_WZTo3LNu_amcatnlo,Skim_TriLep_TTZToLLNuNu -n 10 -r $RUN --userflags Run3Mu
elif [[ $RUN == "Run3" ]]; then
    SKNano.py -a ClosDiLepTrigs -i Skim_TriLep_WZTo3LNu_powheg,Skim_TriLep_TTZ_M50,TTZ_M4to50 -n 10 -r $RUN --userflags Run1E2Mu
    SKNano.py -a ClosDiLepTrigs -i Skim_TriLep_WZTo3LNu_powheg,Skim_TriLep_TTZ_M50,TTZ_M4to50 -n 10 -r $RUN --userflags Run3Mu
else
    echo "Unknown run: $RUN"
    exit 1
fi
