#ifndef PromptTreeProducer_h
#define PromptTreeProducer_h

#include "TriLeptonBase.h"

class PromptTreeProducer : public TriLeptonBase {
public:
    PromptTreeProducer();
    ~PromptTreeProducer();

    void initializeAnalyzer() override;
    void executeEvent() override;
    void WriteHist() override;

private:
    // Channel configuration
    enum class Channel { NONE, SR1E2Mu, SR3Mu };
    Channel channel;

    // Systematic configuration
    std::unique_ptr<SystematicHelper> systHelper;

    // Output trees (one per systematic)
    std::map<TString, TTree*> trees;

    // Tree branches (one set per systematic)
    std::map<TString, float> mass1;
    std::map<TString, float> mass2;
    std::map<TString, float> score_MHc160_MA85_vs_nonprompt;
    std::map<TString, float> score_MHc160_MA85_vs_diboson;
    std::map<TString, float> score_MHc160_MA85_vs_ttZ;
    std::map<TString, float> score_MHc130_MA90_vs_nonprompt;
    std::map<TString, float> score_MHc130_MA90_vs_diboson;
    std::map<TString, float> score_MHc130_MA90_vs_ttZ;
    std::map<TString, float> score_MHc100_MA95_vs_nonprompt;
    std::map<TString, float> score_MHc100_MA95_vs_diboson;
    std::map<TString, float> score_MHc100_MA95_vs_ttZ;
    std::map<TString, int> fold;
    std::map<TString, float> weight;

    // Helper structures
    struct RecoObjects {
        RVec<Muon> vetoMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> vetoElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> jets;
        RVec<Jet> bjets;
        Particle METv;
    };

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float muonRecoSF;
        float muonIDSF;
        float eleRecoSF;
        float eleIDSF;
        float trigSF;
        float pileupIDSF;
        float btagSF;
        float WZNjetsSF;
    };

    // Core methods
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons,
                             const RVec<Electron>& rawElectrons,
                             const RVec<Jet>& rawJets,
                             const RVec<GenJet>& genJets,
                             const TString& syst);

    Channel selectEvent(Event& ev, const RecoObjects& recoObjects,
                       const RVec<Gen>& truth, const TString& syst);

    WeightInfo getWeights(Channel channel, const Event& event,
                         const RecoObjects& recoObjects,
                         const RVec<GenJet>& genJets,
                         const TString& syst);

    void fillTree(Channel channel, const RecoObjects& recoObjects,
                  const WeightInfo& weights, const TString& syst,
                  const Particle& centralMETv);

    void processWeightOnlySystematics(Channel channel, const Event& event,
                                     const RecoObjects& recoObjects,
                                     const RVec<GenJet>& genJets,
                                     const Particle& centralMETv);

    // Helper methods
    Particle makePair(const RVec<Muon>& muons);
    std::pair<Particle, Particle> makePairs(const RVec<Muon>& muons);
    void initTreeContents();
    int calculateFold(const Particle& centralMETv);
    void evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                   const RVec<Jet>& jets, const RVec<Jet>& bjets,
                   const Particle& METv, const TString& syst);
};

#endif
