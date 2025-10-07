#ifndef EvtTreeProducer_h
#define EvtTreeProducer_h

#include "AnalyzerCore.h"

class EvtTreeProducer: public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();
    void WriteHist();

    EvtTreeProducer();
    ~EvtTreeProducer();

    bool Run1E2Mu, Run3Mu;
    IDContainer *MuonIDs, *ElectronIDs;
    RVec<TString> DblMuTriggers, EMuTriggers;

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

    enum class Channel {
        NONE,
        SR1E2MU,
        SR3MU
    };

    RecoObjects defineObjects(const Event& ev, 
                              const RVec<Muon>& rawMuons, 
                              const RVec<Electron>& rawElectrons, 
                              const RVec<Jet>& rawJets, 
                              const RVec<GenJet>& genJets, 
                              const TString& syst);
    Channel selectEvent(const Event& ev, 
                        const RecoObjects& recoObjects, 
                        const TString& syst);

private:
    TTree *newtree;

    // event weights
    float genWeight;
    float puWeight;
    float prefireWeight;

    // muons
    unsigned int nMuons;
    float MuonPtColl[3];
    float MuonEtaColl[3];
    float MuonPhiColl[3];
    float MuonMassColl[3];
    int MuonChargeColl[3];
    bool MuonLabelColl[3];

    // electrons
    unsigned int nElectrons;
    float ElectronPtColl[1];
    float ElectronEtaColl[1];
    float ElectronPhiColl[1];
    float ElectronMassColl[1];
    int ElectronChargeColl[1];
    bool ElectronLabelColl[1];

    // jets
    unsigned int nJets;
    unsigned int nBJets;
    float JetPtColl[20];
    float JetEtaColl[20];
    float JetPhiColl[20];
    float JetMassColl[20];
    float JetChargeColl[20];
    float JetBtagScoreColl[20];
    bool JetIsBtaggedColl[20];
    bool JetLabelColl[20];

    // METv
    float METvPt, METvPhi;
};

#endif
