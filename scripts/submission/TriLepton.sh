#!/bin/bash
RUN=$1
CHANNEL=$2
MODE=$3

## Update DATASTREAM
DATASTREAM=""
MC_LIST=""
if [[ $RUN == "Run2"]]; then
    if [[ $CHANNEL == "Run1E2Mu" || $CHANNEL == "Run2E1Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    MC_LIST="SampleLists/Run2NanoV9"
elif [[ $RUN == "Run3"]]; then
    if [[ $CHANNEL == "Run1E2Mu" || $CHANNEL == "Run2E1Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon,Skim_TriLep_Muon,Skim_TriLep_Muon0,Skim_TriLep_Muon1"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    MC_LIST="SampleLists/Run3NanoV13"
else
    echo "Unknown run: $RUN"
    exit 1
fi

## SR Mode
if [[ $MODE == "SR"]]; then
    SKNano.py -a PromptAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a MatrixAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${CHANNEL} --python
    SKNano.py -a PromptAnalyzer -i ${MC_LIST}/TriLepton.txt -n 20 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 6000
    SKNano.py -a PromptAnalyzer -i ${MC_LIST}/SignalSamples.txt -n 20 -r ${RUN} --userflags ${CHANNEL},RunSyst,RunTheoryUnc --python --memory 6000
elif [[ $MODE == "CR" ]]; then
    SKNano.py -a PromptAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${CHANNEL},RunCR,NoTreeMode --python
    SKNano.py -a MatrixAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${CHANNEL},RunCR,NoTreeMode --python
    SKNano.py -a PromptAnalyzer -i ${MC_LIST}/TriLepton.txt -n 20 -r ${RUN} --userflags ${CHANNEL},RunSyst,RunCR,NoTreeMode --python --memory 6000
    SKNano.py -a PromptAnalyzer -i ${MC_LIST}/SignalSamples.txt -n 20 -r ${RUN} --userflags ${CHANNEL},RunSyst,RunTheoryUnc,RunCR,NoTreeMode --python --memory 6000
elif [[ $MODE == "WZ" ]]; then
    SKNano.py -a PromptAnalyzer -i Skim_TriLep_WZTo3LNu_powheg -n 20 -r ${RUN} --userflags ${CHANNEL},RunSyst --python --memory 6000
    SKNano.py -a PromptAnalyzer -i Skim_TriLep_WZTo3LNu_powheg -n 20 -r ${RUN} --userflags ${CHANNEL},RunNoWZSF,RunSyst --python --memory 6000
else
    echo "Unknown mode: $MODE"
    exit 1
fi
