#ifndef DiLeptonBase_h
#define DiLeptonBase_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class DiLeptonBase: public AnalyzerCore {
public:
    DiLeptonBase();
    ~DiLeptonBase();

    // For di-lepton control regions
    bool RunDiMu, RunEMu;
    bool Run1E2Mu, Run3Mu;
    bool RunNoJetVeto;
    bool RunNoHEMVeto;
    // For fake measurement
    bool MeasFakeMu8, MeasFakeMu17;
    bool MeasFakeEl8, MeasFakeEl12, MeasFakeEl23;
    bool RunSyst;
    
    // IDs
    //RVec<TString> MuonIDs, ElectronIDs;
    IDContainer *MuonIDs, *ElectronIDs;
    // Trigger
    RVec<TString> DblMuTriggers, EMuTriggers;
    TString isoMuTriggerName;
    float triggerSafePtCut;

    //unique_ptr<SystematicHelper> systHelper;
    void initializeAnalyzer();
    virtual void executeEvent();
};

#endif
