#ifndef SignalKinematics_h
#define SignalKinematics_h

#include "TriLeptonBase.h"

class SignalKinematics : public TriLeptonBase {
public:
    SignalKinematics();
    virtual ~SignalKinematics();

    void initializeAnalyzer();
    void executeEvent();

public:
    enum class Channel {
        NONE,
        SR1E2MU,
        SR3MU
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::SR1E2MU) return "SR1E2Mu";
        if (ch == Channel::SR3MU) return "SR3Mu";
        return "NONE";
    }

    struct RecoObjects {
        RVec<Muon> vetoMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> vetoElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> jets;
        RVec<Jet> bjets;
        Particle METv;
    };

private:
    Channel channel;

    RecoObjects defineObjects(Event& ev,
                              const RVec<Muon>& rawMuons,
                              const RVec<Electron>& rawElectrons,
                              const RVec<Jet>& rawJets);

    Channel selectEvent(Event& ev, const RecoObjects& recoObjects);

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float totWeight;
    };

    WeightInfo getWeights(Event& ev, const RecoObjects& recoObjects);

    void fillObjects(Channel ch, const RecoObjects& recoObjects, const WeightInfo& weights, const RVec<Gen>& truth);
    void fillGenLeptonOrigin(const TString& channelStr, Channel ch, const RecoObjects& recoObjects, const RVec<Gen>& truth, float w);
    void fillGenBJetMatching(const TString& channelStr, const RecoObjects& recoObjects, const RVec<Gen>& truth, float w);
    void fillGenDistributions(const RVec<Gen>& truth, float w);
    void fillInclusiveGen(const RVec<Gen>& truth, const RVec<GenJet>& genJets, float w);
    void fillPairSelectionStudy(const TString& channelStr, const RecoObjects& recoObjects, const RVec<Gen>& truth, float w);
};

#endif
