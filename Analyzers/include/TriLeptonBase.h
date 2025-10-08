#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"
#include <TRandom3.h>
#include <map>

// Save ROOT's ClassDef macro before including torch headers
#pragma push_macro("ClassDef")
#pragma push_macro("ClassDefOverride")
#pragma push_macro("ClassDefInline")
#pragma push_macro("ClassDefNV")
#pragma push_macro("ClassImp")

// Undefine ROOT macros temporarily
#undef ClassDef
#undef ClassDefOverride
#undef ClassDefInline
#undef ClassDefNV
#undef ClassImp

// Include torch headers
#include <torch/script.h>
#include <torch/torch.h>

// Restore ROOT macros
#pragma pop_macro("ClassImp")
#pragma pop_macro("ClassDefNV")
#pragma pop_macro("ClassDefInline")
#pragma pop_macro("ClassDefOverride")
#pragma pop_macro("ClassDef")

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
    // ParticleNet models (one per signal point)
    std::map<TString, torch::jit::script::Module> GraphNetModels;
    bool ModelPerFoldWarningIssued;

    // ParticleNet helper functions
    void loadGraphNetModels();
    int calculateFold(const Particle& centralMETv, int nJets);
    torch::Tensor constructNodeFeatures(
        const RVec<Muon*>& muons,
        const RVec<Electron*>& electrons,
        const RVec<Jet*>& jets,
        const RVec<Jet*>& bjets,
        const Particle& METv
    );
    torch::Tensor constructKNNGraph(const torch::Tensor& x, int k = 4);
    torch::Tensor getGraphFeatures(TString era);
    std::map<TString, std::vector<float>> evalGraphNetScores(
        const RVec<Muon*>& muons,
        const RVec<Electron*>& electrons,
        const RVec<Jet*>& jets,
        const RVec<Jet*>& bjets,
        const Particle& METv,
        TString era
    );
};

#endif
