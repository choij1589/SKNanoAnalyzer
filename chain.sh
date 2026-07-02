#!/bin/bash
# DiLepton Validation Region
#./scripts/submission/DiLepton.sh Run2 
#./scripts/submission/DiLepton.sh Run3

# Fake rate measurement
#ERAs=("2016preVFP" "2016postVFP" "2017" "2018"
#      "2022" "2022EE" "2023" "2023BPix")
#OBJECTs=("electron" "muon")

#for ERA in ${ERAs[@]}; do
#    for OBJECT in ${OBJECTs[@]}; do
#        ./scripts/submission/MeasFakeRate.sh $ERA $OBJECT
#    done
#done

#for ERA in ${ERAs[@]}; do
#    for OBJECT in ${OBJECTs[@]}; do
#        ./scripts/convertLeptonFakeRateToJson.py -e $ERA -l $OBJECT
#    done
#done

#SKNano.py -a ClosFakeRate -i Skim_TriLep_TTLL_powheg -n 50 -r Run2,Run3 --userflags Run1E2Mu
#SKNano.py -a ClosFakeRate -i Skim_TriLep_TTLL_powheg -n 50 -r Run2,Run3 --userflags Run3Mu

# Signal Kinematics
#SKNano.py -a SignalKinematics -i SampleLists/Run2NanoV9/SignalSamples.txt  -r Run2 -n 10 --userflags Run1E2Mu
#SKNano.py -a SignalKinematics -i SampleLists/Run2NanoV9/SignalSamples.txt  -r Run2 -n 10 --userflags Run3Mu
#SKNano.py -a SignalKinematics -i SampleLists/Run3NanoV13/SignalSamples.txt  -r Run3 -n 10 --userflags Run1E2Mu
#SKNano.py -a SignalKinematics -i SampleLists/Run3NanoV13/SignalSamples.txt  -r Run3 -n 10 --userflags Run3Mu

# Control regions
#./scripts/submission/TriLepton.sh Run2 Run1E2Mu CR
#./scripts/submission/TriLepton.sh Run2 Run3Mu CR

## In this step, Run3 WZTo3LNu_powheg nominal sample run with invalid WZ NjSF json config
## WZ sample will be run with RunNoWZSF mode
#./scripts/submission/TriLepton.sh Run3 Run1E2Mu CR
#./scripts/submission/TriLepton.sh Run3 Run3Mu CR

## After re-measure WZNjSF for Run3, re-run WZ sample
#SKNano.py -a PromptAnalyzer -i Skim_TriLep_WZTo3LNu_powheg -n 30 -r Run3 --userflags Run1E2Mu,RunSyst,RunCR,NoTreeMode --python --memory 4000
#SKNano.py -a PromptAnalyzer -i Skim_TriLep_WZTo3LNu_powheg -n 30 -r Run3 --userflags Run3Mu,RunSyst,RunCR,NoTreeMode --python --memory 4000

# Signal Region
#./scripts/submission/TriLepton.sh Run2 Run1E2Mu SR
./scripts/submission/TriLepton.sh Run2 Run3Mu SR
./scripts/submission/TriLepton.sh Run2 Run2E1Mu SR
./scripts/submission/TriLepton.sh Run3 Run1E2Mu SR
./scripts/submission/TriLepton.sh Run3 Run3Mu SR
./scripts/submission/TriLepton.sh Run3 Run2E1Mu SR

#./scripts/submission/TriLepton_PN.sh Run2 Run1E2Mu MHc115
#./scripts/submission/TriLepton_PN.sh Run2 Run3Mu MHc115
#./scripts/submission/TriLepton_PN.sh Run2 Run2E1Mu MHc115
#./scripts/submission/TriLepton_PN.sh Run3 Run1E2Mu MHc115
#./scripts/submission/TriLepton_PN.sh Run3 Run3Mu MHc115
#./scripts/submission/TriLepton_PN.sh Run3 Run2E1Mu MHc115
