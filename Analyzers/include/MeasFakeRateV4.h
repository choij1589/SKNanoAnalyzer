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

    // TriggerInfo struct: stores trigger info including event-level fired status
    struct TriggerInfo {
        TString name;           // HLT path: "HLT_Mu8_TrkIsoVVL"
        TString prefix;         // Histogram prefix: "Mu8"
        float trigSafePtCut;    // pT threshold (10, 15, 20, 25 GeV)
        bool fired;             // Event-level: did this trigger fire?
    };

    struct RecoObjects {
        RVec<Muon> vetoMuons;                // 10 GeV, loose ID (for event veto)
        RVec<Muon> looseMuons;               // 10 GeV, loose ID (for measurement)
        RVec<Muon> tightMuons;               // 10 GeV, tight ID (for measurement)
        RVec<Electron> vetoElectrons;        // 10 GeV, loose ID (for event veto)
        RVec<Electron> looseElectrons;       // 15 GeV, loose ID (for measurement)
        RVec<Electron> tightElectrons;       // 15 GeV, tight ID (for measurement)
        RVec<int> looseMuonJetFlavours;      // mother jet hadron flavour for loose muons
        RVec<int> tightMuonJetFlavours;      // mother jet hadron flavour for tight muons
        RVec<int> looseElectronJetFlavours;  // mother jet hadron flavour for loose electrons
        RVec<int> tightElectronJetFlavours;  // mother jet hadron flavour for tight electrons
        RVec<Jet> tightJets;
        RVec<Jet> tightJets_noPUID;  // Jets before PUID for SF calculation (Run 2)
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
    bool MeasFakeMu, MeasFakeEl;
    bool RunSyst;
    bool RunNoHEMVeto;

    // Analysis configuration
    LeptonType leptonType;

    // Triggers for the chosen lepton type
    RVec<TriggerInfo> triggers;

    // Loosest pT cut across all triggers (for initial event selection)
    float lowestPtCut;

    // Binning
    RVec<float> ptcorr_bins;
    RVec<float> abseta_bins;

    // IDs
    IDContainer *MuonIDs, *ElectronIDs;

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
    template<typename T>
    int getMotherJetFlavour(const T& lep, const RVec<Jet>& allJets);
    TString getFlavourSubdir(int flavour);

    // Cutflow functionality
    enum class CutStage {
        Initial = 0,
        NoiseFilter = 1,
        VetoMap = 2,
        AnyTrigger = 3,
        LeptonSelection = 4,
        JetRequirements = 5,
        AwayJetRequirements = 6,
        ZMassWindow = 7,
        Final = 8
    };

    void fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst);
};

#endif
