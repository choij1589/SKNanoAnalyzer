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
        SR3MU
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::SR3MU) return "SR3Mu";
        return "NONE";
    }

    struct RecoObjects {
        RVec<Muon> tightMuons;
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
        float muonRecoSF;
        float muonIDSF;
        float trigSF;
        float btagSF;
        float totWeight;
    };

    WeightInfo getWeights(Event& ev, const RecoObjects& recoObjects);

    void fillObjects(Channel ch, const RecoObjects& recoObjects, const WeightInfo& weights, const RVec<Gen>& truth);
};

#endif
