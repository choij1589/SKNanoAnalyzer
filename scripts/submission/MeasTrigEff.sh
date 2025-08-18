#!/bin/bash
ERA=$1

SKNano.py -a MeasTrigEff -i SingleElectron,EGamma,EGamma0,EGamma1 -n 10 -r $ERA --userflags MeasMuLegs
SKNano.py -a MeasTrigEff -i SingleMuon,Muon,Muon0,Muon1 -n 10 -r $ERA --userflags MeasElLegs
SKNano.py -a MeasTrigEff -i DoubleMuon,Muon,Muon0,Muon1 -n 10 -r $ERA --userflags MeasDblMuPairwise
SKNano.py -a MeasTrigEff -i MuonEG -n 10 -r $ERA --userflags MeasEMuPairwise
SKNano.py -a MeasTrigEff -i TTLL_powheg,DYJets,DYJets10to50 -n 30 -r $ERA --userflags MeasMuLegs
SKNano.py -a MeasTrigEff -i TTLL_powheg,DYJets,DYJets10to50 -n 30 -r $ERA --userflags MeasElLegs
SKNano.py -a MeasTrigEff -i TTLL_powheg,DYJets,DYJets10to50 -n 30 -r $ERA --userflags MeasDblMuPairwise
SKNano.py -a MeasTrigEff -i TTLL_powheg,DYJets,DYJets10to50 -n 30 -r $ERA --userflags MeasEMuPairwise
