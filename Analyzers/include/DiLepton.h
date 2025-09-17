#ifndef DiLepton_h
#define DiLepton_h

#include "DiLeptonBase.h"
#include "MyCorrection.h"
#include "SystematicHelper.h"

class DiLepton: public DiLeptonBase {
public:
    DiLepton();
    virtual ~DiLepton();

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
        RVec<Muon> vetoMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> vetoElectrons;
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
        float muonRecoSF;
        float muonIDSF;
        float eleRecoSF;
        float eleIDSF;
        float trigSF;
        float pileupIDSF;
        float btagSF;
    };

private:

    // channel configuration
    Channel channel;
    bool run_syst;

    // SystematicHelper
    std::unique_ptr<SystematicHelper> systHelper;
    void processWeightOnlySystematics(const Channel& channel, const Event& event, const RecoObjects& recoObjects, const RVec<Gen>& genParts);

    // private methods
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& syst);
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons, 
                             const RVec<Electron>& rawElectrons, 
                             const RVec<Jet>& rawJets,
                             const RVec<GenJet>& genJets, 
                             const TString& syst = "Central");
    WeightInfo getWeights(const Channel& channel,
                          const Event& event,
                          const RecoObjects& recoObjects, 
                          const RVec<Gen>& genParts,
                          const TString& syst = "Central");
    
    void fillObjects(const Channel& channel,
                     const RecoObjects& recoObjects, 
                     const WeightInfo& weights, 
                     const TString& syst = "Central");

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
    
    void fillCutflow(CutStage stage, const Channel& channel, float weight, const TString& syst);
};

#endif
