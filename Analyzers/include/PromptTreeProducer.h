#ifndef PromptTreeProducer_h
#define PromptTreeProducer_h

#include "TriLeptonBase.h"
#include "SystematicHelper.h"

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

    // Theory uncertainty configuration
    RVec<TString> theorySystematics;
    bool hasTheoryWeights;

    // Output trees (one per systematic)
    std::map<TString, TTree*> trees;

    // Tree branches (one set per systematic)
    std::map<TString, float> mass1;
    std::map<TString, float> mass2;
    std::map<TString, float> MT1;
    std::map<TString, float> MT2;
    // Multi-class GraphNet scores: [systematic][masspoint][class] = score
    // masspoints: MHc160_MA85, MHc130_MA90, MHc100_MA95
    // classes: signal, nonprompt, diboson, ttZ
    std::map<TString, std::map<TString, std::map<TString, float>>> ParticleNetScores;
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
    std::tuple<Particle, Particle, float, float> makePairs(const RVec<Muon>& muons, const Particle& METv);
    void initTreeContents();
    void evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                   const RVec<Jet>& jets, const RVec<Jet>& bjets,
                   const Particle& METv, const TString& syst);

    // Theory uncertainty methods
    void initializeTheorySystematics();
    void processTheorySystematics(Channel channel, const RecoObjects& recoObjects,
                                  const WeightInfo& weights, const Particle& centralMETv);
    float getTheoryWeight(const TString& systName);
};

#endif
