#ifndef ClosFakeRate_h
#define ClosFakeRate_h

#include "TriLeptonBase.h"
#include "MyCorrection.h"
#include "SystematicHelper.h"

class ClosFakeRate : public TriLeptonBase {
public:
    ClosFakeRate();
    virtual ~ClosFakeRate();

    void initializeAnalyzer();
    void executeEvent();

public:
    enum class Channel {
        NONE,
        SR1E2MU,     // Signal region 1E2Mu
        SB1E2MU,     // Sideband 1E2Mu  
        SR3MU,       // Signal region 3Mu
        SB3MU        // Sideband 3Mu
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::SR1E2MU) return "SR1E2Mu";
        if (ch == Channel::SB1E2MU) return "SB1E2Mu";
        if (ch == Channel::SR3MU) return "SR3Mu";
        if (ch == Channel::SB3MU) return "SB3Mu";
        return "NONE";
    }

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

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float triggerLumiWeight;
        float fakeWeight;
        float totalWeight;
    };

private:
    // Helper methods
    RecoObjects defineObjects(const Event& ev, 
                              const RVec<Muon>& rawMuons, 
                              const RVec<Electron>& rawElectrons, 
                              const RVec<Jet>& rawJets, 
                              const TString& syst);
    
    Channel selectEvent(const Event& ev, 
                        const RVec<Gen>& truth,
                        const RecoObjects& objects, 
                        const TString& syst);

    WeightInfo getWeights(const Channel selectedChannel,
                          const Event& ev,
                          const RecoObjects& objects,
                          const RVec<Gen>& genParts,
                          const TString& syst);

    void fillObjects(const Channel selectedChannel,
                     const RecoObjects& objects,
                     const WeightInfo& weights,
                     const TString& syst);

    // Helper for 3Mu charge configuration
    std::tuple<Muon, Muon, Muon> configureChargeOf(const RVec<Muon>& muons);
};

#endif