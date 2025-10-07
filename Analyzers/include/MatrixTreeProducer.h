#ifndef MatrixTreeProducer_h
#define MatrixTreeProducer_h

#include "TriLeptonBase.h"

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
    float score_MHc160_MA85_vs_nonprompt;
    float score_MHc160_MA85_vs_diboson;
    float score_MHc160_MA85_vs_ttZ;
    float score_MHc130_MA90_vs_nonprompt;
    float score_MHc130_MA90_vs_diboson;
    float score_MHc130_MA90_vs_ttZ;
    float score_MHc100_MA95_vs_nonprompt;
    float score_MHc100_MA95_vs_diboson;
    float score_MHc100_MA95_vs_ttZ;
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
    std::pair<Particle, Particle> makePairs(const RVec<Muon>& muons);
    void initTreeContents();
    void evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                   const RVec<Jet>& jets, const RVec<Jet>& bjets,
                   const Particle& METv);
};

#endif
