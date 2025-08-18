#ifndef SKNanoLoader_h
#define SKNanoLoader_h

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TROOT.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"
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
    TString Campaign;
    //bool IsFastSim;
    int DataYear;
    TString DataEra;
    TString DataPeriod;
    int Run;
    float xsec, sumW, sumSign;
    RVec<TString> Userflags;

    virtual void Init();
    virtual void SetMaxLeafSize();
    virtual void Loop();
    virtual void executeEvent(){};

    virtual void SetEra(TString era) {
        DataEra=era;
        DataYear=TString(era(0,4)).Atoi();
        //2016, 2017, 2018 are Run2
        if(DataYear == 2016 or DataYear == 2017 or DataYear == 2018) Run = 2;
        else if(DataYear == 2022 or DataYear == 2023 or DataYear == 2024 or DataYear == 2025) Run = 3;
    }
    virtual void SetPeriod(TString period) { DataPeriod=period; }
    virtual void SetCampaign(TString campaign) { Campaign=campaign; }

    virtual TString GetEra() const { return DataEra; }
    virtual TString GetCampaign() const { return Campaign; }
    virtual int GetYear() const { return DataYear; }

    TChain *fChain=nullptr;
    // ****** WARNING ******
    // Do not mix RVec and Array while declaring variables
    // It will cause a memory allocation issue - not be spotted while compiling
    // ****** WARNING ******
    // Weights
    Float_t genWeight;
    Float_t LHEWeight_originalXWGTUP;
    Float_t Generator_weight;
    Int_t nLHEPdfWeight;   // 325300 - 325402
    Int_t nLHEScaleWeight; // 9
    Int_t nPSWeight;       // 4
    RVec<Float_t> LHEPdfWeight;
    RVec<Float_t> LHEScaleWeight;
    RVec<Float_t> PSWeight;

    // PDFs
    Int_t Generator_id1;
    Int_t Generator_id2;
    Float_t Generator_x1;
    Float_t Generator_x2;
    Float_t Generator_xpdf1;
    Float_t Generator_xpdf2;
    Float_t Generator_scalePDF; 

    // LHE
    Float_t LHE_HT, LHE_HTIncoming;
    Float_t LHE_Vpt;
    Float_t LHE_AlphaS;
    UChar_t LHE_Njets, LHE_Nb, LHE_Nc, LHE_Nuds, LHE_Nglu;
    UChar_t LHE_NpLO, LHE_NpNLO;

    // LHEPart
    Int_t nLHEPart;
    UInt_t nLHEPart_RunII;
    RVec<Float_t> LHEPart_pt;
    RVec<Float_t> LHEPart_eta;
    RVec<Float_t> LHEPart_phi;
    RVec<Float_t> LHEPart_mass;
    RVec<Float_t> LHEPart_incomingpz;
    RVec<Int_t> LHEPart_pdgId;
    RVec<Int_t> LHEPart_status;
    RVec<Int_t> LHEPart_spin;
    
    // GenPart
    Int_t nGenPart;
    UInt_t nGenPart_RunII;
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

    // GenJet
    Int_t nGenJet;
    UInt_t nGenJet_RunII;
    RVec<Float_t> GenJet_eta;
    RVec<UChar_t> GenJet_hadronFlavour;
    RVec<Float_t> GenJet_mass;
    RVec<Float_t> GenJet_phi;
    RVec<Float_t> GenJet_pt;
    //Run3
    RVec<Short_t> GenJet_partonFlavour;
    //Run2
    RVec<Int_t> GenJet_partonFlavour_RunII;

    // GenJetAK8
    Int_t nGenJetAK8;
    UInt_t nGenJetAK8_RunII;
    RVec<Float_t> GenJetAK8_eta;
    RVec<UChar_t> GenJetAK8_hadronFlavour;
    RVec<Float_t> GenJetAK8_mass;
    RVec<Float_t> GenJetAK8_phi;
    RVec<Float_t> GenJetAK8_pt;
    //Run3
    RVec<Short_t> GenJetAK8_partonFlavour;
    //Run2
    RVec<Int_t> GenJetAK8_partonFlavour_RunII;

    // GenMet
    Float_t GenMET_pt;
    Float_t GenMET_phi;

    // GenDressedLepton
    Int_t nGenDressedLepton;
    UInt_t nGenDressedLepton_RunII;
    RVec<Float_t> GenDressedLepton_pt;
    RVec<Float_t> GenDressedLepton_eta;
    RVec<Float_t> GenDressedLepton_phi;
    RVec<Float_t> GenDressedLepton_mass;
    RVec<Int_t> GenDressedLepton_pdgId;
    RVec<Bool_t> GenDressedLepton_hasTauAnc;

    // GenIsolatedPhotons
    Int_t nGenIsolatedPhoton;
    UInt_t nGenIsolatedPhoton_RunII;
    RVec<Float_t> GenIsolatedPhoton_pt;
    RVec<Float_t> GenIsolatedPhoton_eta;
    RVec<Float_t> GenIsolatedPhoton_phi;
    RVec<Float_t> GenIsolatedPhoton_mass;
    
    // GenVisTau
    Int_t nGenVisTau;
    UInt_t nGenVisTau_RunII;
    RVec<Float_t> GenVisTau_pt;
    RVec<Float_t> GenVisTau_eta;
    RVec<Float_t> GenVisTau_phi;
    RVec<Float_t> GenVisTau_mass;
    RVec<Int_t> GenVisTau_charge;
    RVec<Int_t> GenVisTau_genPartIdxMother;
    RVec<Int_t> GenVisTau_status;

    // GenVtx -> Need Update

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

    // Muon----------------------------
    Int_t nMuon;
    UInt_t nMuon_RunII;
    RVec<Int_t> Muon_charge;
    RVec<Float_t> Muon_dxy;
    RVec<Float_t> Muon_dxyErr;
    RVec<Float_t> Muon_dxybs;
    RVec<Float_t> Muon_dz;
    RVec<Float_t> Muon_dzErr;
    RVec<Float_t> Muon_eta;
    RVec<UChar_t> Muon_highPtId;
    RVec<Float_t> Muon_ip3d;
    RVec<Int_t> Muon_nTrackerLayers;
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
    RVec<UChar_t> Muon_genPartFlav;
    // Run3
    RVec<UChar_t> Muon_mvaMuID_WP;
    RVec<Short_t> Muon_genPartIdx;
    RVec<Short_t> Muon_jetIdx;
    // Run2
    RVec<UChar_t> Muon_mvaId; //this is in fact wp
    RVec<Int_t> Muon_genPartIdx_RunII;
    RVec<Int_t> Muon_jetIdx_RunII;

    //Electron----------------------------
    Int_t nElectron;
    UInt_t nElectron_RunII;
    RVec<Int_t> Electron_charge;
    RVec<Bool_t> Electron_convVeto;
    RVec<Bool_t> Electron_cutBased_HEEP;
    RVec<Float_t> Electron_scEta;
    RVec<Float_t> Electron_deltaEtaInSC;
    RVec<Float_t> Electron_deltaEtaInSeed;
    RVec<Float_t> Electron_deltaPhiInSC;
    RVec<Float_t> Electron_deltaPhiInSeed;
    RVec<Float_t> Electron_ecalPFClusterIso;
    RVec<Float_t> Electron_hcalPFClusterIso;
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
    RVec<UChar_t> Electron_genPartFlav;
    //Run3
    RVec<UChar_t> Electron_cutBased;
    RVec<Short_t> Electron_fsrPhotonIdx;
    RVec<Short_t> Electron_genPartIdx;
    RVec<Short_t> Electron_jetIdx;
    RVec<Float_t> Electron_mvaHZZIso;
    RVec<Float_t> Electron_mvaIso;
    RVec<Bool_t> Electron_mvaIso_WP80;
    RVec<Bool_t> Electron_mvaIso_WP90;
    RVec<Bool_t> Electron_mvaIso_WPL;
    RVec<Float_t> Electron_mvaNoIso;
    RVec<Bool_t> Electron_mvaNoIso_WP80;
    RVec<Bool_t> Electron_mvaNoIso_WP90;
    RVec<Bool_t> Electron_mvaNoIso_WPL;
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
    RVec<Bool_t> Electron_mvaFall17V2Iso_WPL;
    RVec<Float_t> Electron_mvaFall17V2noIso;
    RVec<Bool_t> Electron_mvaFall17V2noIso_WP80;
    RVec<Bool_t> Electron_mvaFall17V2noIso_WP90;
    RVec<Bool_t> Electron_mvaFall17V2noIso_WPL;
    RVec<Float_t> Electron_dEsigmaUp;
    RVec<Float_t> Electron_dEsigmaDown;

    //Photon----------------------------
    Int_t nPhoton;
    UInt_t nPhoton_RunII;
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
    Int_t nJet;
    UInt_t nJet_RunII;
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
    RVec<UInt_t> Jet_chMultiplicity;
    RVec<UInt_t> Jet_neMultiplicity;
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
    Int_t nTau;
    UInt_t nTau_RunII;
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
    Int_t nFatJet;
    UInt_t nFatJet_RunII;
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

    // MET
    Float_t MET_pt;
    Float_t MET_phi;
    Float_t PuppiMET_pt, PuppiMET_ptUnclusteredUp, PuppiMET_ptUnclusteredDown;
    Float_t PuppiMET_phi, PuppiMET_phiUnclusteredUp, PuppiMET_phiUnclusteredDown;

    // rho
    Float_t fixedGridRhoFastjetAll;
    // Flag
    Bool_t Flag_METFilters; // What is this?
    Bool_t Flag_goodVertices;
    Bool_t Flag_globalSuperTightHalo2016Filter;
    Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter;
    Bool_t Flag_HBHENoiseFilter;
    Bool_t Flag_HBHENoiseIsoFilter;
    Bool_t Flag_BadPFMuonFilter;
    Bool_t Flag_BadPFMuonDzFilter;
    Bool_t Flag_hfNoisyHitsFilter;
    Bool_t Flag_ecalBadCalibFilter;
    Bool_t Flag_eeBadScFilter;
    Int_t RunNumber;
    Int_t LumiBlock;
    Int_t EventNumber;

    // TrigObj----------------------------
    Int_t nTrigObj;
    UInt_t nTrigObj_RunII;
    RVec<Float_t> TrigObj_pt;
    RVec<Float_t> TrigObj_eta;
    RVec<Float_t> TrigObj_phi;
    RVec<Int_t> TrigObj_id_RunII;
    RVec<UShort_t> TrigObj_id;
    RVec<Int_t> TrigObj_filterBits;
    std::map<TString, pair<Bool_t*,float>> TriggerMap;
};

#endif
