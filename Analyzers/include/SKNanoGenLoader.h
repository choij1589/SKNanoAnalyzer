#ifndef SKNanoGenLoader_h
#define SKNanoGenLoader_h

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

class SKNanoGenLoader {
public:
    SKNanoGenLoader();
    virtual ~SKNanoGenLoader();

    virtual void SetTreeName(TString tname) { fChain = new TChain(tname); }
    virtual int AddFile(TString filename) { return fChain->Add(filename, -1); }

    long MaxEvent, NSkipEvent;
    int LogEvery;
    TString MCSample;

    TString Campaign;
    float xsec, sumW, sumSign;
    RVec<TString> Userflags;
    virtual void SetCampaign(TString campaign) { Campaign = campaign; };

    virtual void Init();
    virtual void SetMaxLeafSize();
    virtual void Loop();

    virtual void executeEvent() {};

    virtual TString GetCampaign() const {return Campaign;}

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
    Float_t Generator_weight;
    Float_t Generator_scalePDF;
    Float_t Generator_binvar;

    // LHE
    // LHE scale variation weights (w_var / w_nominal); [0] is MUF="0.5" MUR="0.5"; [1] is MUF="1.0" MUR="0.5"; [2] is MUF="2.0" MUR="0.5"; [3] is MUF="0.5" MUR="1.0"; [4] is MUF="2.0" MUR="1.0"; [5] is MUF="0.5" MUR="2.0"; [6] is MUF="1.0" MUR="2.0"; [7] is MUF="2.0" MUR="2.0"

    Float_t LHEWeight_originalXWGTUP;
    static constexpr int nLHEPdfWeight = 110; // 325300 - 325402
    static constexpr int nLHEScaleWeight = 9;
    Float_t LHEPdfWeight[nLHEPdfWeight];
    Float_t LHEScaleWeight[nLHEScaleWeight];
    Float_t LHE_HT, LHE_HTIncoming;
    Float_t LHE_Vpt;
    Float_t LHE_AlphaS;
    UChar_t LHE_Njets, LHE_Nb, LHE_Nc, LHE_Nuds, LHE_Nglu;
    UChar_t LHE_NpLO, LHE_NpNLO;

    UInt_t* nLHEPart = new UInt_t;
    RVec<Float_t> LHEPart_pt;
    RVec<Float_t> LHEPart_eta;
    RVec<Float_t> LHEPart_phi;
    RVec<Float_t> LHEPart_mass;
    RVec<Float_t> LHEPart_incomingpz;
    RVec<Int_t> LHEPart_pdgId;
    RVec<Int_t> LHEPart_status;
    RVec<Int_t> LHEPart_spin;

    // GenPart
    // GenPart statusFlags bits:
    // 0=isPrompt, 1=isDecayedLeptonHadron, 2=isTauDecayProduct,
    // 3=isPromptTauDecayProduct, 4=isDirectTauDecayProduct,
    // 5=isDirectPromptTauDecayProduct, 6=isDirectHadronDecayProduct,
    // 7=isHardProcess, 8=fromHardProcess, 9=isHardProcessTauDecayProduct,
    // 10=isDirectHardProcessTauDecayProduct, 11=fromHardProcessBeforeFSR,
    // 12=isFirstCopy, 13=isLastCopy, 14=isLastCopyBeforeFSR
    // For the mass, quarks(except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411) are not stored. For these particles, you can lookup the value from PDG.
    UInt_t* nGenPart = new UInt_t;
    RVec<Float_t> GenPart_pt;
    RVec<Float_t> GenPart_eta;
    RVec<Float_t> GenPart_phi;
    RVec<Float_t> GenPart_mass;
    RVec<Int_t> GenPart_pdgId;
    RVec<Int_t> GenPart_status;
    RVec<Int_t> GenPart_statusFlags;
    RVec<Int_t> GenPart_genPartIdxMother;

    // GenJet
    UInt_t* nGenJet = new UInt_t;
    RVec<Float_t> GenJet_pt;
    RVec<Float_t> GenJet_eta;
    RVec<Float_t> GenJet_phi;
    RVec<Float_t> GenJet_mass;
    RVec<Int_t> GenJet_partonFlavour;
    RVec<UChar_t> GenJet_hadronFlavour;

    UInt_t* nGenJetAK8 = new UInt_t;
    RVec<Float_t> GenJetAK8_pt;
    RVec<Float_t> GenJetAK8_eta;
    RVec<Float_t> GenJetAK8_phi;
    RVec<Float_t> GenJetAK8_mass;
    RVec<Int_t> GenJetAK8_partonFlavour;
    RVec<UChar_t> GenJetAK8_hadronFlavour;

    // GenMET
    Float_t GenMET_pt;
    Float_t GenMET_phi;
    Float_t MET_fiducialGenPt;
    Float_t MET_fiducialGenPhi;

    // GenDressedLeptons
    UInt_t* nGenDressedLepton = new UInt_t;
    RVec<Float_t> GenDressedLepton_pt;
    RVec<Float_t> GenDressedLepton_eta;
    RVec<Float_t> GenDressedLepton_phi;
    RVec<Float_t> GenDressedLepton_mass;
    RVec<Int_t> GenDressedLepton_pdgId;
    RVec<Bool_t> GenDressedLepton_hasTauAnc;

    // GenIsolatedPhotons
    UInt_t* nGenIsolatedPhoton = new UInt_t;
    RVec<Float_t> GenIsolatedPhoton_pt;
    RVec<Float_t> GenIsolatedPhoton_eta;
    RVec<Float_t> GenIsolatedPhoton_phi;
    RVec<Float_t> GenIsolatedPhoton_mass;

    // VisTau
    // VisTau status bits 
    // 0=OneProng0PiZero, 1=OneProng1PiZero, 2=OneProng2PiZero
    // 10=ThreeProng0PiZero, 11=ThreeProng1PiZero, 15=Other
    UInt_t* nGenVisTau = new UInt_t;
    RVec<Float_t> GenVisTau_pt;
    RVec<Float_t> GenVisTau_eta;
    RVec<Float_t> GenVisTau_phi;
    RVec<Float_t> GenVisTau_mass;
    RVec<Int_t> GenVisTau_charge;
    RVec<Int_t> GenVisTau_genPartIdxMother;
    RVec<Int_t> GenVisTau_status;

    // GenVtx -> Need Update
};

#endif