#ifndef MeasFakeRateV4_h
#define MeasFakeRateV4_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "MyCorrection.h"

class MeasFakeRateV4: public AnalyzerCore {
public:
    MeasFakeRateV4();
    virtual ~MeasFakeRateV4();

    void initializeAnalyzer();
    void executeEvent();

    enum class Channel {
        NONE,
        INCLUSIVE,
        QCDENRICHED,
        WENRICHED,
        ZENRICHED
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::INCLUSIVE) return "Inclusive";
        if (ch == Channel::QCDENRICHED) return "QCDEnriched";
        if (ch == Channel::WENRICHED) return "WEnriched";
        if (ch == Channel::ZENRICHED) return "ZEnriched";
        return "NONE";
    }

    enum class LeptonType {
        NONE,
        MUON,
        ELECTRON
    };

    struct RecoObjects {
        RVec<Muon> looseMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> looseElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> tightJets;
        //RVec<Jet> tightJets_vetoLep;
        RVec<Jet> bjets;
        RVec<GenJet> genJets;
        Particle METv;
    };

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float topPtWeight;
        float muonRecoSF;
        float eleRecoSF;
        float btagSF;
        float pileupIDSF;
    };

private:
    // Configuration flags
    bool MeasFakeMu8, MeasFakeMu17;
    bool MeasFakeEl8, MeasFakeEl12, MeasFakeEl23;
    bool MeasFakeMu, MeasFakeEl;
    bool RunSyst;

    // Analysis configuration
    LeptonType leptonType;
    
    // Binning
    RVec<float> ptcorr_bins;
    RVec<float> abseta_bins;
    
    // IDs
    IDContainer *MuonIDs, *ElectronIDs;
    
    // Trigger
    TString isoSglLepTrig;
    float trigSafePtCut;

    // SystematicHelper
    std::unique_ptr<SystematicHelper> systHelper;

    // Core analysis methods
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst);
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons, 
                             const RVec<Electron>& rawElectrons, 
                             const RVec<Jet>& rawJets,
                             const RVec<GenJet>& genJets,
                             const TString& ID,
                             const TString& syst = "Central");
    WeightInfo getWeights(const Channel& channel,
                          const TString& ID,
                          const Event& event,
                          const RecoObjects& recoObjects, 
                          const RVec<Gen>& genParts,
                          const TString& syst = "Central");
    
    void fillObjects(const Channel& channel,
                     const TString& ID,
                     const RecoObjects& recoObjects, 
                     const WeightInfo& weights, 
                     const TString& syst = "Central");

    // Helper methods
    TString getBinPrefix(const double ptcorr, const double abseta);
    float getJetPtCut(const TString& selection);

    // Cutflow functionality
    enum class CutStage {
        Initial = 0,
        NoiseFilter = 1,
        VetoMap = 2,
        Trigger = 3,
        LeptonSelection = 4,
        JetRequirements = 5,
        AwayJetRequirements = 6,
        ZMassWindow = 7,
        Final = 8
    };
    
    void fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst);
};

#endif
