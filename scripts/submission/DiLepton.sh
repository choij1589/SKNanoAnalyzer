#!/bin/bash
ERA=$1
CHANNEL=("DiMu" "EMu")

# If ERA == "201*", use Run2NanoV9
# If ERA == "202*", use Run3NanoV13


for ch in "${CHANNEL[@]}"; do
    echo "Submitting jobs for channel: $ch"

    if [[ $ch == "EMu" ]]; then
        DATASTREAM="MuonEG"
    else
        DATASTREAM="DoubleMuon,Muon,Muon0,Muon1"
    fi

    if [[ $ERA == "Run2" ]]; then
        SKNano.py -a DiLepton -i $DATASTREAM -n 20 -r $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i DYJets,WJets,TTLL_powheg,TTLJ_powheg -n 50 -r $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i SampleLists/Run2NanoV9/DiLepton.txt -n 20 -r $ERA --userflags Run${ch},RunSyst
    elif [[ $ERA == "Run3" ]]; then
        SKNano.py -a DiLepton -i $DATASTREAM -n 20 -r $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i DYJets,WJets,TTLL_powheg,TTLJ_powheg -n 50 -r $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i SampleLists/Run3NanoV13/DiLepton.txt -n 20 -r $ERA --userflags Run${ch},RunSyst
    elif [[ $ERA == "201"* ]]; then
        SKNano.py -a DiLepton -i $DATASTREAM -n 20 -e $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i DYJets,WJets,TTLL_powheg,TTLJ_powheg -n 50 -e $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i SampleLists/Run2NanoV9/DiLepton.txt -n 20 -e $ERA --userflags Run${ch},RunSyst
    elif [[ $ERA == "202"* ]]; then
        SKNano.py -a DiLepton -i $DATASTREAM -n 20 -e $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i DYJets,WJets,TTLL_powheg,TTLJ_powheg -n 50 -e $ERA --userflags Run${ch},RunSyst
        SKNano.py -a DiLepton -i SampleLists/Run3NanoV13/DiLepton.txt -n 20 -e $ERA --userflags Run${ch},RunSyst
    fi
done
