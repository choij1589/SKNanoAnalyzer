#ifndef MeasFakeRate_h
#define MeasFakeRate_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "MyCorrection.h"

class MeasFakeRate: public AnalyzerCore {
public:
    MeasFakeRate();
    virtual ~MeasFakeRate();

    void initializeAnalyzer();
    void executeEvent();

public:
    enum class Channel {
        NONE,
        INCLUSIVE,
        ZENRICHED
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::INCLUSIVE) return "Inclusive";
        if (ch == Channel::ZENRICHED) return "ZEnriched";
        return "NONE";
    }

    enum class LeptonType {
        MUON,
        ELECTRON
    };

    struct RecoObjects {
        RVec<Muon> looseMuons;
        RVec<Muon> tightMuons;
        RVec<Muon> vetoMuons;
        RVec<Electron> looseElectrons;
        RVec<Electron> tightElectrons;
        RVec<Electron> vetoElectrons;
        RVec<Jet> tightJets;
        RVec<Jet> tightJets_vetoLep;
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
    };

private:
    // Configuration flags
    bool MeasFakeMu8, MeasFakeMu17;
    bool MeasFakeEl8, MeasFakeEl12, MeasFakeEl23;
    bool MeasFakeMu, MeasFakeEl;
    bool RunSyst;

    // Analysis configuration
    LeptonType leptonType;
    TString currentID;
    
    // Binning
    RVec<double> ptcorr_bins;
    RVec<double> abseta_bins;
    
    // IDs
    IDContainer *MuonIDs, *ElectronIDs;
    
    // Trigger
    TString isoSglLepTrig;
    float trigSafePtCut;

    // SystematicHelper
    std::unique_ptr<SystematicHelper> systHelper;
    void processWeightOnlySystematics(const Channel& channel, const TString& ID, const Event& event, const RecoObjects& recoObjects, const RVec<Gen>& genParts);

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
    TString findBin(const double ptcorr, const double abseta);
    double getJetPtCut(const TString& selection);
};

#endif