#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"
#include <TRandom3.h>
#include <map>

// Torch headers removed - using Python PyTorch for ParticleNet inference

class TriLeptonBase: public AnalyzerCore {
public:
    TriLeptonBase();
    ~TriLeptonBase();

    // For tri-lepton regions
    bool Run1E2Mu, Run3Mu, Run2E1Mu;
    bool RunSyst;
    bool RunTheoryUnc;
    bool RunNoVetoMap;
    bool RunNoWZSF;

    // IDs
    IDContainer *MuonIDs, *ElectronIDs;
    // Trigger
    RVec<TString> DblMuTriggers, EMuTriggers;

    void initializeAnalyzer();
    virtual void executeEvent();

    float GetFakeWeight(const RVec<Muon> &muons, const RVec<Electron> &electrons, const TString syst_key="Central");
    RVec<Electron> GetPTCorrScaledElectrons(const RVec<Electron> &electros);
    RVec<Muon> GetPTCorrScaledMuons(const RVec<Muon> &muons);
protected:
    // Fold calculation (still useful for validation)
    int calculateFold(const Particle& centralMETv, int nJets);
};

#endif
