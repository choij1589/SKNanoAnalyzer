#!/bin/bash
ERA=$1
OBJECT=$2

DATASTREAM=""
SAMPLELIST=""
QCD=""
##### DATASTREAM and QCD
if [[ $OBJECT == "muon" ]]; then
    QCD="QCD_MuEnriched"
    if [[ $ERA == "201"* ]]; then
        DATASTREAM="DoubleMuon"
        SAMPLELIST="SampleLists/Run2NanoV9"
    elif [[ $ERA == "2022" ]]; then
        DATASTREAM="DoubleMuon,Muon"
        SAMPLELIST="SampleLists/Run3NanoV13"
    elif [[ $ERA == "2022EE" ]]; then
        DATASTREAM="Muon"
        SAMPLELIST="SampleLists/Run3NanoV13"
    elif [[ $ERA == "2023"* ]]; then
        DATASTREAM="Muon0,Muon1"
        SAMPLELIST="SampleLists/Run3NanoV13"
    else 
        echo "Wrong era $ERA"
        exit 1
    fi
    USERFLAGs=("MeasFakeMu8" "MeasFakeMu17")
    for USERFLAG in ${USERFLAGs[@]}; do
        SKNano.py -a MeasFakeRateV2 -i $DATASTREAM -n 10 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i DYJets,WJets,TTLJ_powheg,TTLL_powheg -n 30 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i $SAMPLELIST/DiLepton.txt -n 10 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i $SAMPLELIST/$QCD.txt -n 10 -e $ERA --userflags $USERFLAG,RunSyst
    done
elif [[ $OBJECT == "electron" ]]; then
    QCD="QCD_EMEnriched"
    if [[ $ERA == "2016"* ]]; then
        DATASTREAM="DoubleEG"
        SAMPLELIST="SampleLists/Run2NanoV9"
    elif [[ $ERA == "2017" ]]; then
        DATASTREAM="SingleElectron"
        SAMPLELIST="SampleLists/Run2NanoV9"
    elif [[ $ERA == "2018" ]]; then
        DATASTREAM="EGamma"
        SAMPLELIST="SampleLists/Run2NanoV9"
    elif [[ $ERA == "2022"* ]]; then
        DATASTREAM="EGamma"
        SAMPLELIST="SampleLists/Run3NanoV13"
    elif [[ $ERA == "2023"* ]]; then
        DATASTREAM="EGamma0,EGamma1"
        SAMPLELIST="SampleLists/Run3NanoV13"
    else
        echo "Wrong era $ERA"
        exit 1
    fi
    USERFLAGs=("MeasFakeEl8" "MeasFakeEl12" "MeasFakeEl23")
    for USERFLAG in ${USERFLAGs[@]}; do
        SKNano.py -a MeasFakeRateV2 -i $DATASTREAM -n 10 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i DYJets,WJets,TTLJ_powheg,TTLL_powheg -n 30 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i $SAMPLELIST/DiLepton.txt -n 10 -e $ERA --userflags $USERFLAG,RunSyst
        SKNano.py -a MeasFakeRateV2 -i $SAMPLELIST/$QCD.txt -n 10 -e $ERA --userflags $USERFLAG,RunSyst 
    done
else
    echo "Wrong object $OBJECT"
    exit 1
fi
