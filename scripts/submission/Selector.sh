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
    SKNano.py -a PromptSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a MatrixSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a PromptSelector -i SampleLists/Run2NanoV9/TriLepton.txt -n 10 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 12000
    SKNano.py -a PromptSelector -i SampleLists/Run2NanoV9/SignalSamples.txt -n 10 -r ${RUN} --userflags ${CHANNEL},RunSyst,RunTheoryUnc --python --memory 12000
elif [[ $RUN == "Run3" ]]; then
    if [[ $CHANNEL == "Run1E2Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon,Skim_TriLep_Muon,Skim_TriLep_Muon0,Skim_TriLep_Muon1"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    SKNano.py -a PromptSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a MatrixSelector -i $DATASTREAM -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a PromptSelector -i SampleLists/Run3NanoV13/TriLepton.txt -n 10 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 12000
fi
