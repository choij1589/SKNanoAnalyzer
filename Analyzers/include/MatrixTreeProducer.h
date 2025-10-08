#ifndef MatrixTreeProducer_h
#define MatrixTreeProducer_h

#include "TriLeptonBase.h"
#include "SystematicHelper.h"

class MatrixTreeProducer : public TriLeptonBase {
public:
    MatrixTreeProducer();
    ~MatrixTreeProducer();

    void initializeAnalyzer() override;
    void executeEvent() override;
    void WriteHist() override;

private:
    // Channel configuration
    enum class Channel { NONE, SR1E2Mu, SR3Mu };
    Channel channel;

    // Output tree
    TTree* newtree;

    // Tree branches
    float mass1;
    float mass2;
    float MT1;
    float MT2;
    // Multi-class GraphNet scores: [masspoint][class] = score
    // masspoints: MHc160_MA85, MHc130_MA90, MHc100_MA95
    // classes: signal, nonprompt, diboson, ttZ
    std::map<TString, std::map<TString, float>> ParticleNetScores;
    int fold;
    float weight;

    // Helper structures
    struct RecoObjects {
        RVec<Muon> vetoMuons;
        RVec<Muon> looseMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> vetoElectrons;
        RVec<Electron> looseElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> jets;
        RVec<Jet> bjets;
        Particle METv;
    };

    // Core methods
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons,
                             const RVec<Electron>& rawElectrons,
                             const RVec<Jet>& rawJets);

    Channel selectEvent(Event& ev, const RecoObjects& recoObjects);

    void fillTree(Channel channel, const RecoObjects& recoObjects);

    // Helper methods
    Particle makePair(const RVec<Muon>& muons);
    std::tuple<Particle, Particle, float, float> makePairs(const RVec<Muon>& muons, const Particle& METv);
    void initTreeContents();
    void evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                   const RVec<Jet>& jets, const RVec<Jet>& bjets,
                   const Particle& METv);
};

#endif
