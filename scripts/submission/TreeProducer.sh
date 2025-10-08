#!/bin/bash
ERA=$1
CHANNEL=$2

SAMPLELIST=""
#if [[ $RUN == "Run2" ]]; then
SAMPLELIST="SampleLists/Run2NanoV9"
#elif [[ $RUN == "Run3" ]]; then
#    SAMPLELIST="SampleLists/Run3NanoV13"
#else
#    echo "Unknown run: $RUN"
#    exit 1
#fi

DATASTREAM=""
if [[ $CHANNEL == "Run1E2Mu" ]]; then
    DATASTREAM="Skim_TriLep_MuonEG"
elif [[ $CHANNEL == "Run3Mu" ]]; then
    DATASTREAM="Skim_TriLep_DoubleMuon,Skim_TriLep_Muon,Skim_TriLep_Muon0,Skim_TriLep_Muon1"
else
    echo "Unknown channel: $CHANNEL"
    exit 1
fi

SKNano.py -a PromptTreeProducer -i $DATASTREAM -n 10 -e $ERA --userflags $CHANNEL
SKNano.py -a PromptTreeProducer -i $SAMPLELIST/TriLepton.txt -n 10 -e $ERA --userflags $CHANNEL,RunSyst
SKNano.py -a MatrixTreeProducer -i $DATASTREAM -n 10 -e $ERA --userflags $CHANNEL
SKNano.py -a PromptTreeProducer -i $SAMPLELIST/SignalSamples.txt -n 1 -e $ERA --userflags $CHANNEL,RunSyst
