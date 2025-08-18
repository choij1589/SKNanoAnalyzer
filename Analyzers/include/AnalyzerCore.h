#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include <map>
#include <unordered_map>
#include <string>
#include <deque>

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObjString.h"
#include "TMath.h"

#include "SKNanoLoader.h"
#include "Event.h"
#include "Particle.h"
#include "Lepton.h"
#include "Gen.h"
#include "LHE.h"
#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "FatJet.h"
#include "Tau.h"
#include "Photon.h"
#include "Gen.h"
#include "Jet.h"
#include "GenJet.h"
#include "GenDressedLepton.h"
#include "GenIsolatedPhoton.h"
#include "GenVisTau.h"
#include "TrigObj.h"

#include "LHAPDFHandler.h"
#include "PDFReweight.h"
#include "MyCorrection.h"
#include "JetTaggingParameter.h"
#include "PhysicalConstants.h"

class IDContainer {
public:
    IDContainer() {}
    IDContainer(const TString &tight, const TString &loose):
        j_tight(tight), j_loose(loose) {}

    TString GetID(const TString &wp) const {
        if (wp == "tight") return j_tight;
        else if (wp == "loose") return j_loose;
        else throw runtime_error("Invalid WP: " + wp);
    }

private:
    TString j_tight, j_loose;
};

class AnalyzerCore: public SKNanoLoader {
public:
    AnalyzerCore();
    ~AnalyzerCore();

    virtual void initializeAnalyzer() {};
    virtual void executeEvent() {};

    inline bool HasFlag(const TString &flag) { return std::find(Userflags.begin(), Userflags.end(), flag) != Userflags.end(); }

    inline static bool PtComparing(const Particle& p1, const Particle& p2) { return p1.Pt() > p2.Pt();}
    inline static bool PtComparingPtr(const Particle* p1, const Particle* p2) { return p1->Pt() > p2->Pt();}


    //MetFilter
    bool PassNoiseFilter(const RVec<Jet> &AllJets, const Event &ev, const Event::MET_Type MET_type=Event::MET_Type::PUPPI);
    // PDF reweight
    PDFReweight *pdfReweight;
    float GetPDFWeight(LHAPDF::PDF *pdf_);
    float GetPDFReweight();
    float GetPDFReweight(int member);    
    // Correction
    MyCorrection *myCorr;
    //unique_ptr<CorrectionSet> csetMuon;
    //unique_ptr<CorrectionSet> csetElectron;;

    // MC weights
    float MCweight(bool usesign = true, bool norm_1invpb = true) const;

    // Get objects
    Event GetEvent();
    RVec<Muon> GetAllMuons();
    RVec<Muon> GetMuons(const TString ID, const float ptmin, const float fetamax);
    RVec<Electron> GetAllElectrons();
    RVec<Jet> GetAllJets();
    RVec<Gen> GetAllGens();
    RVec<LHE> GetAllLHEs();
    RVec<Jet> GetJets(const TString id, const float ptmin, const float fetamax);
    RVec<Electron> GetElectrons(const TString id, const float ptmin, const float fetamax, bool vetoHEM = false);
    RVec<Tau> GetAllTaus();
    RVec<FatJet> GetAllFatJets();
    RVec<GenJet> GetAllGenJets();
    RVec<GenDressedLepton> GetAllGenDressedLeptons();
    RVec<GenIsolatedPhoton> GetAllGenIsolatedPhotons();
    RVec<GenVisTau> GetAllGenVisTaus();
    RVec<Photon> GetAllPhotons();
    RVec<Photon> GetPhotons(TString id, double ptmin, double fetamax);
    RVec<TrigObj> GetAllTrigObjs();

    // Select objects
    RVec<Muon> SelectMuons(const RVec<Muon> &muons, TString ID, const float ptmin, const float absetamax) const;
    RVec<Muon> SelectMuons(const RVec<Muon> &muons, Muon::MuonID ID, const float ptmin, const float absetamax) const;
    RVec<Jet> SelectJets(const RVec<Jet> &jets, const TString id, const float ptmin, const float fetamax) const;
    RVec<Jet> SelectJets(const RVec<Jet> &jets, const Jet::JetID, const float ptmin, const float fetamax) const;
    RVec<Jet> JetsVetoLeptonInside(const RVec<Jet> &jets, const RVec<Electron> &electrons, const RVec<Muon> &muons, const float dR = 0.3) const;
    RVec<Electron> SelectElectrons(const RVec<Electron> &electrons, const TString id, const float ptmin, const float absetamax, bool vetoHEM = false) const;
    RVec<Electron> SelectElectrons(const RVec<Electron> &electrons, const Electron::ElectronID ID, const float ptmin, const float absetamax, bool vetoHEM = false) const;
    RVec<Tau> SelectTaus(const RVec<Tau> &taus, const TString ID, const float ptmin, const float absetamax) const;
    // Functions
    float GetScaleVariation(const MyCorrection::variation &muF_syst, const MyCorrection::variation &muR_syst);
    float GetPSWeight(const MyCorrection::variation &ISR_syst, const MyCorrection::variation &FSR_syst);
    inline float GetBTaggingWP(const JetTagging::JetFlavTagger &tagger, const JetTagging::JetFlavTaggerWP &wp) { return myCorr->GetBTaggingWP(tagger, wp); }
    inline pair<float, float> GetCTaggingWP(const JetTagging::JetFlavTagger &tagger, const JetTagging::JetFlavTaggerWP &wp) { return myCorr->GetCTaggingWP(tagger, wp); }
    inline float GetBTaggingWP(){ return myCorr->GetBTaggingWP(); }
    inline pair<float, float> GetCTaggingWP(){ return myCorr->GetCTaggingWP(); }
    float GetHT(const RVec<Jet> &jets);
    bool IsHEMElectron(const Electron& electron) const;

    // Gen Matching
    void PrintGen(const RVec<Gen> &gens);
    static RVec<int> TrackGenSelfHistory(const Gen& me, const RVec<Gen>& gens);
    static Gen GetGenMatchedLepton(const Lepton& lep, const RVec<Gen>& gens);
    static Gen GetGenMatchedMuon(const Muon& muon, const RVec<Gen>& gens);
    static Gen GetGenMatchedPhoton(const Lepton& lep, const RVec<Gen>& gens);
    static bool IsFinalPhotonSt23_Public(const RVec<Gen>& gens);
    bool IsFromHadron(const Gen& me, const RVec<Gen>& gens);
    bool IsSignalPID(const int &pid);
    int GetLeptonType(const Lepton& lep, const RVec<Gen>& gens);
    int GetLeptonType(const Gen& gen, const RVec<Gen>& gens);
    int GetLeptonType_Public(const int& genIdx, const RVec<Gen>& gens);
    int GetGenPhotonType(const Gen& genph, const RVec<Gen>& gens);
    int GetPrElType_InSameSCRange_Public(int genIdx, const RVec<Gen>& gens);

    // Scale and smear
    void METType1Propagation(Particle &MET, RVec<Particle> &original_objects, RVec<Particle> &corrected_objects);
    float GetL1PrefireWeight(MyCorrection::variation syst = MyCorrection::variation::nom);
    unordered_map<int, int> GenJetMatching(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const float &rho, const float dR = 0.2, const float pTJerCut = 3.);
    unordered_map<int, int> deltaRMatching(const RVec<Particle> &objs1, const RVec<Particle> &objs2, const float dR = 0.4);
    RVec<Muon> ScaleMuons(const RVec<Muon> &muons, const TString &syst );
    RVec<Electron> ScaleElectrons(const Event &ev, const RVec<Electron> &electrons, const TString &syst);
    RVec<Electron> SmearElectrons(const RVec<Electron> &electrons, const TString &syst);

    RVec<Jet> SmearJets(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    RVec<Jet> SmearJets(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const TString &syst, const TString &source="total");
    RVec<Jet> ScaleJets(const RVec<Jet> &jets, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    RVec<Jet> ScaleJets(const RVec<Jet> &jets, const TString &syst, const TString &source="total");
    
    // Type-I MET correction with correlated object variations
    Particle ApplyTypeICorrection(const Particle& MET,
                                  const RVec<Jet>& jets,
                                  const RVec<Electron>& electrons,
                                  const RVec<Muon>& muons,
                                  const MyCorrection::variation& unclustered_syst=MyCorrection::variation::nom);
    
    // Histogram Handlers
    TFile* GetOutfile() { return outfile; }
    inline void SetOutfilePath(const TString &outpath) { outfile = new TFile(outpath, "RECREATE"); }
    TH1D* GetHist1D(const string &histname);
    bool PassVetoMap(const Jet &jet, const RVec<Muon> &AllMuons, const TString mapCategory="jetvetomap");
    bool PassVetoMap(const RVec<Jet> &AllJets, const RVec<Muon> &AllMuons, const TString mapCategory="jetvetomap");
    bool PassJetVetoMap(const RVec<Jet> &AllJet, const RVec<Muon> &AllMuon, const TString mapCategory="jetvetomap");
    inline void FillCutFlow(const int &val,const int &maxCutN=10){
        static int storedMaxCutN = maxCutN;
        FillHist("CutFlow", val, 1., storedMaxCutN, 0, storedMaxCutN);
    }
    void FillHist(const TString &histname, float value, float weight, int n_bin, float x_min, float x_max);
    void FillHist(const TString &histname, float value, float weight, int n_bin, float *xbins);
    void FillHist(const TString &histname, float value_x, float value_y, float weight, 
                                          int n_binx, float x_min, float x_max, 
                                          int n_biny, float y_min, float y_max);
    void FillHist(const TString &histname, float value_x, float value_y, float weight,
                                          int n_binx, float *xbins,
                                          int n_biny, float *ybins);
    void FillHist(const TString &histname, float value_x, float value_y, float value_z, float weight,
                                          int n_binx, float x_min, float x_max,
                                          int n_biny, float y_min, float y_max,
                                          int n_binz, float z_min, float z_max);
    void FillHist(const TString &histname, float value_x, float value_y, float value_z, float weight,
                                          int n_binx, float *xbins,
                                          int n_biny, float *ybins,
                                          int n_binz, float *zbins);
    inline void FillHist(const TString &histname, float value, float weight, const RVec<float> &xbins) { FillHist(histname, value, weight, xbins.size()-1, const_cast<float*>(xbins.data())); } 
    inline void FillHist(const TString &histname, float value_x, float value_y, float weight, const RVec<float> &xbins, const RVec<float> &ybins) {FillHist(histname, value_x, value_y, weight, xbins.size() - 1, const_cast<float *>(xbins.data()), ybins.size() - 1, const_cast<float *>(ybins.data())); }
    inline void FillHist(const TString &histname, float value_x, float value_y, float value_z, float weight, const RVec<float> &xbins, const RVec<float> &ybins, const RVec<float> &zbins) {FillHist(histname, value_x, value_y, value_z, weight, xbins.size() - 1, const_cast<float *>(xbins.data()), ybins.size() - 1, const_cast<float *>(ybins.data()), zbins.size() - 1, const_cast<float *>(zbins.data())); }


    TTree* NewTree(const TString &treename, const RVec<TString> &keeps = {}, const RVec<TString> &drops = {});
    TTree* GetTree(const TString &treename);
    inline void SetBranch(const TString &treename, const TString &branchname, float val) { this_floats.push_back(val); SetBranch(treename, branchname, (void*)(&this_floats.back()), branchname + "/F"); };
    inline void SetBranch(const TString &treename, const TString &branchname, double val) { this_floats.push_back(float(val)); SetBranch(treename, branchname, (void*)(&this_floats.back()), branchname + "/F"); };
    inline void SetBranch(const TString &treename, const TString &branchname, int val) { this_ints.push_back(val); SetBranch(treename, branchname, (void*)(&this_ints.back()), branchname + "/I"); };
    inline void SetBranch(const TString &treename, const TString &branchname, bool val) { this_bools.push_back(val); SetBranch(treename, branchname, (void *)(&this_bools.back()), branchname + "/O"); }
    //fill RVec to branch -> Not work do not use
    //template <typename T>
    //inline void SetBranch(const TString &treename, const TString &branchname, std::vector<T> &val) {SetBranch_Vector(treename, branchname, val);};

    void FillTrees(const TString &treename="");
    virtual void WriteHist();

private:
    bool useTH1F = false;
    unordered_map<string, TH1*> histmap1d;
    unordered_map<string, TH2*> histmap2d;
    unordered_map<string, TH3*> histmap3d;
    unordered_map<string, TTree*> treemap;
    unordered_map<TTree*, unordered_map<string, TBranch*>> branchmaps; 
    deque<float> this_floats;
    deque<int> this_ints;
    deque<char> this_bools;
    TFile *outfile;
    void SetBranch(const TString &treename, const TString &branchname, void *address, const TString &leaflist);
    template <typename T>
    void SetBranch_Vector(const TString &treename, const TString &branchname, std::vector<T> &address) {
        //Not work do not use
        try {
            TTree *tree = GetTree(treename);

            unordered_map<string, TBranch *> *this_branchmap = &branchmaps[tree];
            auto it = this_branchmap->find(string(branchname));

            if (it == this_branchmap->end())
            {
                //template <typename T, std::size_t N> TBranch *Branch(const char* name, std::array<T, N> *obj, Int_t bufsize = 32000, Int_t splitlevel = 99)
                auto br = tree->Branch(branchname, &address);
                this_branchmap->insert({string(branchname), br});
            }
            else
            {
                //void TBranch::SetAddress(void *add)
                it -> second->SetAddress(&address);
            }
        } catch (int e) {
            cout << "[AnalyzerCore::SetBranch] Error get tree: " << treename.Data() << endl;
            exit(e);
        }
    }
};

#endif
