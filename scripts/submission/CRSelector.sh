#!/bin/bash
RUN=$1
CHANNEL=$2

if [[ $RUN == "Run2" ]]; then
    if [[ $CHANNEL == "Run1E2Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    SKNano.py -a CRPromptSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a CRMatrixSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a CRPromptSelector -i SampleLists/Run2NanoV9/TriLepton.txt -n 10 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 8000
elif [[ $RUN == "Run3" ]]; then
    if [[ $CHANNEL == "Run1E2Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon,Skim_TriLep_Muon,Skim_TriLep_Muon0,Skim_TriLep_Muon1"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    SKNano.py -a CRPromptSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a CRMatrixSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a CRPromptSelector -i SampleLists/Run3NanoV13/TriLepton.txt -n 10 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 8000
    SKNano.py -a CRPromptSelector -i Skim_TriLep_WZTo3LNu_powheg -n 10 -r ${RUN} --userflags ${CHANNEL},RunNoWZSF,RunSyst --python --memory 8000
fi
