#ifndef SKNanoLoader_h
#define SKNanoLoader_h

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include <nlohmann/json.hpp>
using namespace ROOT::VecOps;

class SKNanoLoader {
public:
    SKNanoLoader();
    virtual ~SKNanoLoader();

    //virtual long GetEntry(long entry);
    virtual void SetTreeName(TString tname) { fChain = new TChain(tname); }
    virtual int AddFile(TString filename) { return fChain->Add(filename, -1); }

    long MaxEvent, NSkipEvent;
    int LogEvery;
    bool IsDATA;
    TString DataStream;
    TString MCSample;
    //bool IsFastSim;
    int DataYear;
    TString DataEra;
    int Run;
    float xsec, sumW, sumSign;
    RVec<TString> Userflags;

    virtual void Init();
    virtual void SetMaxLeafSize();
    virtual void Loop();

    //virtual void beginEvent(){};
    virtual void executeEvent(){};
    //virtual void endEvent(){};

    virtual void SetEra(TString era) {
        DataEra=era;
        DataYear=TString(era(0,4)).Atoi();
        //2016, 2017, 2018 are Run2
        if(DataYear == 2016 or DataYear == 2017 or DataYear == 2018) Run = 2;
        else if(DataYear == 2022 or DataYear == 2023 or DataYear == 2024 or DataYear == 2025) Run = 3;
    }
    
    virtual TString GetEra() const { return DataEra; }
    //virtual TString GetEraShort() const;
    virtual int GetYear() const { return DataYear; }

    TChain *fChain=nullptr;
    ROOT::RDataFrame *df=nullptr;

    // Declaration of leaf types
    // PDFs
    Int_t Generator_id1;
    Int_t Generator_id2;
    Float_t Generator_x1;
    Float_t Generator_x2;
    Float_t Generator_xpdf1;
    Float_t Generator_xpdf2;
    Float_t Generator_scalePDF;
    Float_t Generator_weight;
    // Weights
    Float_t genWeight;
    Float_t LHEWeight_originalXWGTUP;
    static constexpr int nLHEPdfWeight = 110; // 325300 - 325402
    static constexpr int nLHEScaleWeight = 9;
    static constexpr int nPSWeight = 4; 
    Float_t LHEPdfWeight[nLHEPdfWeight];
    Float_t LHEScaleWeight[nLHEScaleWeight];
    Float_t PSWeight[nPSWeight];
    //L1Prefire-------------------------
    Float_t L1PreFiringWeight_Nom;
    Float_t L1PreFiringWeight_Dn;
    Float_t L1PreFiringWeight_Up;
    // Event-------------------------
    Int_t Pileup_nPU;
    Float_t Pileup_nTrueInt;
    Int_t genTtbarId;

    // PV----------------------------
    Float_t PV_chi2;
    Float_t PV_ndof;
    Float_t PV_score;
    Float_t PV_x;
    Float_t PV_y;
    Float_t PV_z;
    // Run3
    UChar_t PV_npvs;
    UChar_t PV_npvsGood;
    // Run2
    Int_t PV_npvs_RunII;
    Int_t PV_npvsGood_RunII;

    //GenPart----------------------------
    Int_t* nGenPart = new Int_t;
    UInt_t* nGenPart_RunII = new UInt_t;
    RVec<Float_t> GenPart_eta;
    RVec<Float_t> GenPart_mass;
    RVec<Int_t> GenPart_pdgId;
    RVec<Float_t> GenPart_phi;
    RVec<Float_t> GenPart_pt;
    RVec<Int_t> GenPart_status;
    //Run3
    RVec<Short_t> GenPart_genPartIdxMother;
    RVec<UShort_t> GenPart_statusFlags;
    //Run2
    RVec<Int_t> GenPart_genPartIdxMother_RunII;
    RVec<Int_t> GenPart_statusFlags_RunII;

    //LHEPart-----------------
    Int_t* nLHEPart = new Int_t;
    UInt_t* nLHEPart_RunII = new UInt_t;
    RVec<Float_t> LHEPart_pt;
    RVec<Float_t> LHEPart_eta;
    RVec<Float_t> LHEPart_phi;
    RVec<Float_t> LHEPart_mass;
    RVec<Float_t> LHEPart_incomingpz;
    RVec<Int_t> LHEPart_pdgId;
    RVec<Int_t> LHEPart_status;
    RVec<Int_t> LHEPart_spin;

    // Muon----------------------------
    Int_t* nMuon = new Int_t;
    UInt_t* nMuon_RunII = new UInt_t;
    RVec<Int_t> Muon_charge;
    RVec<Float_t> Muon_dxy;
    RVec<Float_t> Muon_dxyErr;
    RVec<Float_t> Muon_dxybs;
    RVec<Float_t> Muon_dz;
    RVec<Float_t> Muon_dzErr;
    RVec<Float_t> Muon_eta;
    RVec<UChar_t> Muon_highPtId;
    RVec<Float_t> Muon_ip3d;
    RVec<Bool_t> Muon_isGlobal;
    RVec<Bool_t> Muon_isStandalone;
    RVec<Bool_t> Muon_isTracker;
    RVec<Bool_t> Muon_looseId;
    RVec<Float_t> Muon_mass;
    RVec<Bool_t> Muon_mediumId;
    RVec<Bool_t> Muon_mediumPromptId;
    RVec<UChar_t> Muon_miniIsoId;
    RVec<Float_t> Muon_miniPFRelIso_all;
    RVec<UChar_t> Muon_multiIsoId;
    RVec<Float_t> Muon_mvaLowPt;
    RVec<Float_t> Muon_mvaTTH;
    RVec<UChar_t> Muon_pfIsoId;
    RVec<Float_t> Muon_pfRelIso03_all;
    RVec<Float_t> Muon_pfRelIso04_all;
    RVec<Float_t> Muon_phi;
    RVec<Float_t> Muon_pt;
    RVec<UChar_t> Muon_puppiIsoId;
    RVec<Float_t> Muon_sip3d;
    RVec<Bool_t> Muon_softId;
    RVec<Float_t> Muon_softMva;
    RVec<Bool_t> Muon_softMvaId;
    RVec<Bool_t> Muon_tightId;
    RVec<UChar_t> Muon_tkIsoId;
    RVec<Float_t> Muon_tkRelIso;
    RVec<Bool_t> Muon_triggerIdLoose;
    // Run3
    RVec<UChar_t> Muon_mvaMuID_WP;
    // Run2
    RVec<UChar_t> Muon_mvaId; //this is in fact wp

    //Electron----------------------------
    Int_t* nElectron = new Int_t;
    UInt_t* nElectron_RunII = new UInt_t;
    RVec<Int_t> Electron_charge;
    RVec<Bool_t> Electron_convVeto;
    RVec<Bool_t> Electron_cutBased_HEEP;
    RVec<Float_t> Electron_deltaEtaSC;
    RVec<Float_t> Electron_dr03EcalRecHitSumEt;
    RVec<Float_t> Electron_dr03HcalDepth1TowerSumEt;
    RVec<Float_t> Electron_dr03TkSumPt;
    RVec<Float_t> Electron_dr03TkSumPtHEEP;
    RVec<Float_t> Electron_dxy;
    RVec<Float_t> Electron_dxyErr;
    RVec<Float_t> Electron_dz;
    RVec<Float_t> Electron_dzErr;
    RVec<Float_t> Electron_eInvMinusPInv;
    RVec<Float_t> Electron_energyErr;
    RVec<Float_t> Electron_eta;
    RVec<UChar_t> Electron_genPartFlav;
    RVec<Float_t> Electron_hoe;
    RVec<Float_t> Electron_ip3d;
    RVec<Bool_t> Electron_isPFcand;
    RVec<UChar_t> Electron_jetNDauCharged;
    RVec<Float_t> Electron_jetPtRelv2;
    RVec<Float_t> Electron_jetRelIso;
    RVec<UChar_t> Electron_lostHits;
    RVec<Float_t> Electron_mass;
    RVec<Float_t> Electron_miniPFRelIso_all;
    RVec<Float_t> Electron_miniPFRelIso_chg;
    RVec<Float_t> Electron_mvaTTH;
    RVec<Int_t> Electron_pdgId;
    RVec<Float_t> Electron_pfRelIso03_all;
    RVec<Float_t> Electron_pfRelIso03_chg;
    RVec<Float_t> Electron_phi;
    RVec<Float_t> Electron_pt;
    RVec<Float_t> Electron_r9;
    RVec<Float_t> Electron_scEtOverPt;
    RVec<UChar_t> Electron_seedGain;
    RVec<Float_t> Electron_sieie;
    RVec<Float_t> Electron_sip3d;
    //Run3
    RVec<UChar_t> Electron_cutBased;
    RVec<Short_t> Electron_fsrPhotonIdx;
    RVec<Short_t> Electron_genPartIdx;
    RVec<Short_t> Electron_jetIdx;
    RVec<Float_t> Electron_mvaHZZIso;
    RVec<Float_t> Electron_mvaIso;
    RVec<Bool_t> Electron_mvaIso_WP80;
    RVec<Bool_t> Electron_mvaIso_WP90;
    RVec<Float_t> Electron_mvaNoIso;
    RVec<Bool_t> Electron_mvaNoIso_WP80;
    RVec<Bool_t> Electron_mvaNoIso_WP90;
    RVec<Short_t> Electron_photonIdx;
    RVec<Char_t> Electron_seediEtaOriX;
    RVec<Int_t> Electron_seediPhiOriY;
    RVec<Short_t> Electron_svIdx;
    RVec<UChar_t> Electron_tightCharge;
    //Run2
    RVec<UChar_t> Electron_cleanmask;
    RVec<Int_t> Electron_cutBased_RunII;
    RVec<Int_t> Electron_genPartIdx_RunII;
    RVec<Int_t> Electron_jetIdx_RunII;
    RVec<Float_t> Electron_mvaFall17V2Iso;
    RVec<Bool_t> Electron_mvaFall17V2Iso_WP80;
    RVec<Bool_t> Electron_mvaFall17V2Iso_WP90;
    RVec<Float_t> Electron_mvaFall17V2noIso;
    RVec<Bool_t> Electron_mvaFall17V2noIso_WP80;
    RVec<Bool_t> Electron_mvaFall17V2noIso_WP90;

    //Photon----------------------------
    Int_t* nPhoton = new Int_t;
    UInt_t* nPhoton_RunII = new UInt_t;
    RVec<Float_t> Photon_energyErr;
    RVec<Float_t> Photon_eta;
    RVec<UChar_t> Photon_genPartFlav;
    RVec<Float_t> Photon_hoe;
    RVec<Bool_t> Photon_isScEtaEB;
    RVec<Bool_t> Photon_isScEtaEE;
    RVec<Float_t> Photon_mvaID;
    RVec<Bool_t> Photon_mvaID_WP80;
    RVec<Bool_t> Photon_mvaID_WP90;
    RVec<Float_t> Photon_phi;
    RVec<Bool_t> Photon_pixelSeed;
    RVec<Float_t> Photon_pt;
    RVec<Float_t> Photon_energyRaw;
    RVec<Float_t> Photon_sieie;
    //Run3
    RVec<UChar_t> Photon_cutBased;
    //Run2
    RVec<Int_t> Photon_cutBased_RunII;

    //Jet----------------------------
    Int_t* nJet = new Int_t;
    UInt_t* nJet_RunII = new UInt_t;
    RVec<Float_t> Jet_area;
    RVec<Float_t> Jet_btagDeepFlavB;
    RVec<Float_t> Jet_btagDeepFlavCvB;
    RVec<Float_t> Jet_btagDeepFlavCvL;
    RVec<Float_t> Jet_btagDeepFlavQG;
    RVec<Float_t> Jet_chEmEF;
    RVec<Float_t> Jet_chHEF;
    RVec<Float_t> Jet_eta;
    RVec<Int_t> Jet_hfadjacentEtaStripsSize;
    RVec<Int_t> Jet_hfcentralEtaStripSize;
    RVec<Float_t> Jet_hfsigmaEtaEta;
    RVec<Float_t> Jet_hfsigmaPhiPhi;
    RVec<Float_t> Jet_mass;
    RVec<Float_t> Jet_muEF;
    RVec<Float_t> Jet_muonSubtrFactor;
    RVec<UChar_t> Jet_nConstituents;
    RVec<Float_t> Jet_neEmEF;
    RVec<Float_t> Jet_neHEF;
    RVec<Float_t> Jet_phi;
    RVec<Float_t> Jet_pt;
    RVec<Float_t> Jet_rawFactor;
    //Run3
    RVec<Float_t> Jet_PNetRegPtRawCorr;
    RVec<Float_t> Jet_PNetRegPtRawCorrNeutrino;
    RVec<Float_t> Jet_PNetRegPtRawRes;
    RVec<Float_t> Jet_btagPNetB;
    RVec<Float_t> Jet_btagPNetCvB;
    RVec<Float_t> Jet_btagPNetCvL;
    RVec<Float_t> Jet_btagPNetQvG;
    RVec<Float_t> Jet_btagPNetTauVJet;
    RVec<Float_t> Jet_btagRobustParTAK4B;
    RVec<Float_t> Jet_btagRobustParTAK4CvB;
    RVec<Float_t> Jet_btagRobustParTAK4CvL;
    RVec<Float_t> Jet_btagRobustParTAK4QG;
    RVec<Short_t> Jet_electronIdx1;
    RVec<Short_t> Jet_electronIdx2;
    RVec<Short_t> Jet_genJetIdx;
    RVec<UChar_t> Jet_hadronFlavour;
    RVec<UChar_t> Jet_jetId;
    RVec<Short_t> Jet_muonIdx1;
    RVec<Short_t> Jet_muonIdx2;
    RVec<UChar_t> Jet_nElectrons;
    RVec<UChar_t> Jet_nMuons;
    RVec<UChar_t> Jet_nSVs;
    RVec<Short_t> Jet_partonFlavour;
    RVec<Short_t> Jet_svIdx1;
    RVec<Short_t> Jet_svIdx2;
    //Run2
    RVec<Float_t> Jet_bRegCorr;
    RVec<Float_t> Jet_bRegRes;
    RVec<Float_t> Jet_btagCSVV2;
    //RVec<Float_t> Jet_btagDeepB;
    //RVec<Float_t> Jet_btagDeepCvB;
    //RVec<Float_t> Jet_btagDeepCvL;
    RVec<Float_t> Jet_cRegCorr;
    RVec<Float_t> Jet_cRegRes;
    RVec<Float_t> Jet_chFPV0EF;
    RVec<UChar_t> Jet_cleanmask;
    RVec<Int_t> Jet_electronIdx1_RunII;
    RVec<Int_t> Jet_electronIdx2_RunII;
    RVec<Int_t> Jet_genJetIdx_RunII;
    RVec<Int_t> Jet_hadronFlavour_RunII;
    RVec<Int_t> Jet_jetId_RunII;
    RVec<Int_t> Jet_muonIdx1_RunII;
    RVec<Int_t> Jet_muonIdx2_RunII;
    RVec<Int_t> Jet_nElectrons_RunII;
    RVec<Int_t> Jet_nMuons_RunII;
    RVec<Int_t> Jet_partonFlavour_RunII;
    RVec<Int_t> Jet_puId;
    RVec<Float_t> Jet_puIdDisc;
    RVec<Float_t> Jet_qgl;

    //Tau---------------------------------------
    Int_t* nTau = new Int_t;
    UInt_t* nTau_RunII = new UInt_t;
    RVec<Float_t> Tau_dxy;
    RVec<Float_t> Tau_dz;
    RVec<Float_t> Tau_eta;
    RVec<UChar_t> Tau_genPartFlav;
    RVec<UChar_t> Tau_idDeepTau2017v2p1VSe;
    RVec<UChar_t> Tau_idDeepTau2017v2p1VSjet;
    RVec<UChar_t> Tau_idDeepTau2017v2p1VSmu;
    RVec<Float_t> Tau_mass;
    RVec<Float_t> Tau_phi;
    RVec<Float_t> Tau_pt;
    // Run3
    RVec<Short_t> Tau_charge;
    RVec<UChar_t> Tau_decayMode;
    RVec<Short_t> Tau_genPartIdx;
    RVec<Bool_t> Tau_idDecayModeNewDMs;
    RVec<UChar_t> Tau_idDeepTau2018v2p5VSe;
    RVec<UChar_t> Tau_idDeepTau2018v2p5VSjet;
    RVec<UChar_t> Tau_idDeepTau2018v2p5VSmu;
    // Run2
    RVec<Int_t> Tau_charge_RunII;
    RVec<Int_t> Tau_decayMode_RunII;
    RVec<Int_t> Tau_genPartIdx_RunII;

    // FatJet----------------------------
    Int_t* nFatJet = new Int_t;
    UInt_t* nFatJet_RunII = new UInt_t;
    RVec<Float_t> FatJet_area;
    RVec<Float_t> FatJet_btagDDBvLV2;
    RVec<Float_t> FatJet_btagDDCvBV2;
    RVec<Float_t> FatJet_btagDDCvLV2;
    RVec<Float_t> FatJet_btagDeepB;
    RVec<Float_t> FatJet_btagHbb;
    RVec<Float_t> FatJet_eta;
    RVec<Float_t> FatJet_lsf3;
    RVec<Float_t> FatJet_mass;
    RVec<Float_t> FatJet_msoftdrop;
    RVec<UChar_t> FatJet_nBHadrons;
    RVec<UChar_t> FatJet_nCHadrons;
    RVec<UChar_t> FatJet_nConstituents;
    RVec<Float_t> FatJet_particleNet_QCD;
    RVec<Float_t> FatJet_phi;
    RVec<Float_t> FatJet_pt;
    RVec<Float_t> FatJet_tau1;
    RVec<Float_t> FatJet_tau2;
    RVec<Float_t> FatJet_tau3;
    RVec<Float_t> FatJet_tau4;
    // Run3
    RVec<Short_t> FatJet_genJetAK8Idx;
    RVec<UChar_t> FatJet_jetId;
    RVec<Float_t> FatJet_particleNetWithMass_H4qvsQCD;
    RVec<Float_t> FatJet_particleNetWithMass_HbbvsQCD;
    RVec<Float_t> FatJet_particleNetWithMass_HccvsQCD;
    RVec<Float_t> FatJet_particleNetWithMass_QCD;
    RVec<Float_t> FatJet_particleNetWithMass_TvsQCD;
    RVec<Float_t> FatJet_particleNetWithMass_WvsQCD;
    RVec<Float_t> FatJet_particleNetWithMass_ZvsQCD;
    RVec<Float_t> FatJet_particleNet_QCD0HF;
    RVec<Float_t> FatJet_particleNet_QCD1HF;
    RVec<Float_t> FatJet_particleNet_QCD2HF;
    RVec<Float_t> FatJet_particleNet_XbbVsQCD;
    RVec<Float_t> FatJet_particleNet_XccVsQCD;
    RVec<Float_t> FatJet_particleNet_XggVsQCD;
    RVec<Float_t> FatJet_particleNet_XqqVsQCD;
    RVec<Float_t> FatJet_particleNet_XteVsQCD;
    RVec<Float_t> FatJet_particleNet_XtmVsQCD;
    RVec<Float_t> FatJet_particleNet_XttVsQCD;
    RVec<Float_t> FatJet_particleNet_massCorr;
    RVec<Short_t> FatJet_subJetIdx1;
    RVec<Short_t> FatJet_subJetIdx2;
    // Run2
    // I'll pass deepTag scores
    // RVec<Float_t> FatJet_btagCSVV2;
    // RVec<Float_t> FatJet_deepTagMD_H4qvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_HbbvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_TvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_WvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_ZHbbvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_ZHccvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_ZbbvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_ZvsQCD;
    // RVec<Float_t> FatJet_deepTagMD_bbvsLight;
    // RVec<Float_t> FatJet_deepTagMD_ccvsLight;
    // RVec<Float_t> FatJet_deepTag_H;
    // RVec<Float_t> FatJet_deepTag_QCD;
    // RVec<Float_t> FatJet_deepTag_QCDothers;
    // RVec<Float_t> FatJet_deepTag_TvsQCD;
    // RVec<Float_t> FatJet_deepTag_WvsQCD;
    // RVec<Float_t> FatJet_deepTag_ZvsQCD;
    RVec<Int_t> FatJet_genJetAK8Idx_RunII;
    RVec<Int_t> FatJet_jetId_RunII;
    RVec<Float_t> FatJet_particleNetMD_QCD;
    RVec<Float_t> FatJet_particleNetMD_Xbb;
    RVec<Float_t> FatJet_particleNetMD_Xcc;
    RVec<Float_t> FatJet_particleNetMD_Xqq;
    RVec<Float_t> FatJet_particleNet_H4qvsQCD;
    RVec<Float_t> FatJet_particleNet_HbbvsQCD;
    RVec<Float_t> FatJet_particleNet_HccvsQCD;
    RVec<Float_t> FatJet_particleNet_TvsQCD;
    RVec<Float_t> FatJet_particleNet_WvsQCD;
    RVec<Float_t> FatJet_particleNet_ZvsQCD;
    RVec<Float_t> FatJet_particleNet_mass;
    RVec<Int_t> FatJet_subJetIdx1_RunII;
    RVec<Int_t> FatJet_subJetIdx2_RunII;

    //GenJet----------------------------
    Int_t* nGenJet = new Int_t;
    UInt_t* nGenJet_RunII = new UInt_t;
    RVec<Float_t> GenJet_eta;
    RVec<UChar_t> GenJet_hadronFlavour;
    RVec<Float_t> GenJet_mass;
    RVec<Float_t> GenJet_phi;
    RVec<Float_t> GenJet_pt;
    //Run3
    RVec<Short_t> GenJet_partonFlavour;
    //Run2
    RVec<Int_t> GenJet_partonFlavour_RunII;

    //PuppiMET----------------------------
    Float_t PuppiMET_pt;
    Float_t PuppiMET_phi;
    Float_t PuppiMET_sumEt;
    Float_t PuppiMET_ptJERDown;
    Float_t PuppiMET_phiJERDown;
    Float_t PuppiMET_ptJERUp;
    Float_t PuppiMET_phiJERUp;
    Float_t PuppiMET_ptJESDown;
    Float_t PuppiMET_phiJESDown;
    Float_t PuppiMET_ptJESUp;
    Float_t PuppiMET_phiJESUp;
    Float_t PuppiMET_ptUnclusteredDown;
    Float_t PuppiMET_phiUnclusteredDown;
    Float_t PuppiMET_ptUnclusteredUp;
    Float_t PuppiMET_phiUnclusteredUp;

    //MET----------------------------
    Float_t MET_MetUnclustEnUpDeltaX;
    Float_t MET_MetUnclustEnUpDeltaY;
    Float_t MET_covXX;
    Float_t MET_covXY;
    Float_t MET_covYY;
    Float_t MET_fiducialGenPhi;
    Float_t MET_fiducialGenPt;
    Float_t MET_phi;
    Float_t MET_pt;
    Float_t MET_significance;
    Float_t MET_sumEt;
    Float_t MET_sumPtUnclustered;

    // rho
    Float_t fixedGridRhoFastjetAll;
    // Flag
    Bool_t Flag_METFilters; // What is this?
    Bool_t Flag_goodVertices;
    Bool_t Flag_globalSuperTightHalo2016Filter;
    Bool_t Flag_ECalDeadCellTriggerPrimitiveFilter;
    Bool_t Flag_BadPFMuonFilter;
    Bool_t Flag_BadPFMuonDzFilter;
    Bool_t Flag_hfNoisyHitsFilter;
    Bool_t Flag_ecalBadCalibFilter;
    Bool_t Flag_eeBadScFilter;
    //Bool_t Flag_ecalBadCalibFilter;
    UInt_t RunNumber;
    std::map<TString, pair<Bool_t*,float>> TriggerMap;
};

#endif
