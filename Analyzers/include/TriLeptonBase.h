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
    bool Run1E2Mu, Run3Mu;
    bool RunSyst;
    bool RunNoVetoMap;
    bool RunNoWZSF;

    // IDs
    IDContainer *MuonIDs, *ElectronIDs;
    // Trigger
    RVec<TString> DblMuTriggers, EMuTriggers;

    void initializeAnalyzer();
    virtual void executeEvent();

    float GetFakeWeight(const RVec<Muon> &muons, const RVec<Electron> &electrons, const TString syst_key="Central");

protected:
    // ParticleNet C++ implementation removed - use Python ParticleNet instead
    // For Python ParticleNet: See PyAnalyzers/PromptTreeProducer.py
    bool ModelPerFoldWarningIssued;

    // Fold calculation (still useful for validation)
    int calculateFold(const Particle& centralMETv, int nJets);
};

#endif
