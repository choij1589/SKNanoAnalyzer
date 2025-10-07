#ifndef SimpleDiLepton_h
#define SimpleDiLepton_h

#include "DiLeptonBase.h"

class SimpleDiLepton: public DiLeptonBase {
public:
    SimpleDiLepton();
    virtual ~SimpleDiLepton();

    void initializeAnalyzer();
    void executeEvent();

public:
    enum class Channel {
        NONE,
        DIMU,
        EMU
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::DIMU) return "DIMU";
        if (ch == Channel::EMU) return "EMU";
        return "NONE";
    }

    struct RecoObjects {
        RVec<Muon> tightMuons;
        RVec<Electron> tightElectrons;
        RVec<Jet> tightJets;
        RVec<Jet> tightJets_vetoLep;
        RVec<Jet> bjets;
        RVec<GenJet> genJets;
        Particle METv_default;
        Particle METv;
    };

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float topPtWeight;
    };

private:
    // channel configuration
    Channel channel;

    // private methods
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects);
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons,
                             const RVec<Electron>& rawElectrons,
                             const RVec<Jet>& rawJets,
                             const RVec<GenJet>& genJets);
    WeightInfo getWeights(const Channel& channel,
                          const Event& event,
                          const RecoObjects& recoObjects,
                          const RVec<Gen>& genParts);

    void fillObjects(const Channel& channel,
                     const RecoObjects& recoObjects,
                     const WeightInfo& weights);

    // Cutflow functionality
    enum class CutStage {
        Initial = 0,
        NoiseFilter = 1,
        VetoMap = 2,
        LeptonSelection = 3,
        Trigger = 4,
        KinematicCuts = 5,
        JetRequirements = 6,
        BjetRequirements = 7,
        Final = 8
    };

    void fillCutflow(CutStage stage, const Channel& channel, float weight);
};

#endif