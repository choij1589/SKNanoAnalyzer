#ifndef MyCorrection_h
#define MyCorrection_h

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <variant>
using namespace std;

#include "TString.h"
#include "TRandom3.h"
#include "TMath.h"
#include "correction.h"
#include "RoccoR.h"
#include "Jet.h"
#include "JetTaggingParameter.h"
#include "Gen.h"
#include "Muon.h"
#include "Electron.h"
using correction::CorrectionSet;

class MyCorrection {
public:
    enum class variation {
        nom,
        up,
        down
    };

    enum class POG {
        MUO,
        LUM,
        BTV,
        EGM,
        JME
    };

    enum class XYCorrection_MetType {
        Type1PFMET,
        Type1PuppiMET
    };
    MyCorrection();
    MyCorrection(const TString &era, const TString &period, const TString &sample, bool IsData,const string &btagging_eff_file = "btaggingEff.json", const string &ctagging_eff_file = "ctaggingEff.json", const string &btagging_R_file = "btaggingR.json", const string &ctagging_R_file = "ctaggingR.json");
    ~MyCorrection();

    // Muon
    float GetMuonScaleSF(const Muon &muon, const variation syst = variation::nom, const float matched_pt=0) const;
    float GetMuonRECOSF(const Muon &muon, const variation syst = variation::nom) const;
    float GetMuonRECOSF(const RVec<Muon> &muons, const variation syst = variation::nom) const;
    inline float GetMuonISOSF(const TString &Muon_ISO_SF_Key, const Muon &muon, const variation syst = variation::nom, const TString &source = "") { return GetMuonIDSF(Muon_ISO_SF_Key, muon, syst); }
    inline float GetMuonTriggerSF(const TString &Muon_Trigger_SF_Key, const Muon &muon, const variation syst = variation::nom, const TString &source = "") { return GetMuonIDSF(Muon_Trigger_SF_Key, muon, syst); };
    inline float GetMuonISOSF(const TString &Muon_ISO_SF_Key, const RVec<Muon> &muons, const variation syst = variation::nom, const TString &source = "") { return GetMuonIDSF(Muon_ISO_SF_Key, muons, syst); }
    float GetMuonIDSF(const TString &Muon_ID_SF_Key, const Muon &muon, const variation syst = variation::nom) const;
    float GetMuonIDSF(const TString &Muon_ID_SF_Key, const RVec<Muon> &muons, const variation syst = variation::nom) const;
    float GetFakeRate(const Muon muon, const TString &id_key, const TString &syst_key) const;

    // electron
    float GetElectronScaleUnc(const float scEta, const unsigned char seedGain, const unsigned int runNumber, const float r9, const float pt, const variation syst = variation::nom) const;
    float GetElectronSmearUnc(const Electron &electron, const variation syst = variation::nom, const unsigned int seed=999) const;
    float GetElectronRECOSF(const float abseta, const float pt, const float phi, const variation syst = variation::nom) const;
    float GetElectronRECOSF(const RVec<Electron> &electrons, const variation syst = variation::nom) const;
    float GetElectronIDSF(const TString &Electron_ID_SF_Key, const float abseta, const float pt, const float phi, const variation syst = variation::nom) const;
    float GetElectronIDSF(const TString &Electron_ID_SF_Key, const RVec<Electron> &electrons, const variation syst = variation::nom) const;
    float GetFakeRate(const Electron &electron, const TString &id_key, const TString &fakerate_key) const;
    
    // photon

    // Trigger
    // Single lepton trigger from POG
    float GetMuonTriggerEff(const TString &Muon_Trigger_Eff_Key, const float abseta, const float pt, const bool isData, const variation syst = variation::nom) const;
    float GetMuonTriggerSF(const TString &Muon_Trigger_Eff_Key, const RVec<Muon> &muons, const variation syst = variation::nom) const;
    float GetElectronTriggerEff(const TString &Electron_ID_SF_Key, const float eta, const float pt, const float phi, const bool isDATA, const variation syst = variation::nom) const;
    inline float GetElectronTriggerDataEff(const TString &Electron_ID_SF_Key, const float eta, const float pt, const float phi, const variation syst = variation::nom) {
        return GetElectronTriggerEff(Electron_ID_SF_Key, eta, pt, phi, true, syst);
    };
    inline float GetElectronTriggerMCEff(const TString &Electron_ID_SF_Key, const float eta, const float pt, const float phi, const variation syst = variation::nom) {
        return GetElectronTriggerEff(Electron_ID_SF_Key, eta, pt, phi, false, syst);
    };
    float GetElectronTriggerSF(const TString &Electron_Trigger_SF_Key, const float eta, const float pt, const float phi, const variation syst = variation::nom) const;

    // double lepton triggers
    float GetTriggerEff(const Muon &muon, const TString &trigger_leg_key, const bool isData, const variation syst=variation::nom) const;
    float GetTriggerEff(const Electron &electron, const TString &trigger_leg_key, const bool isData, const variation syst=variation::nom) const;
    float GetDblMuTriggerEff(const RVec<Muon> &muons, const bool isDATA, const variation syst=variation::nom) const;
    float GetDblMuTriggerSF(const RVec<Muon> &muons, const variation syst=variation::nom) const;
    float GetEMuTriggerEff(const RVec<Electron> &electrons, const RVec<Muon> &muons, const bool isDATA, const variation syst=variation::nom) const;
    float GetEMuTriggerSF(const RVec<Electron> &electrons, const RVec<Muon> &muons, const variation syst=variation::nom) const;
    float GetPairwiseFilterEff(const TString &filter_name, const bool isData) const;


    // PUWeights
    float GetPUWeight(const float nTrueInt, const variation syst = variation::nom, const TString &source = "") const;

    // tagging param
    void SetTaggingParam(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp);

    // btaging
    float GetBTaggingWP() const;
    float GetBTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const;
    float GetBTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst = variation::nom);
    float GetBTaggingReweightMethod1a(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method, const variation syst=variation::nom, const TString &source="central");
    float GetBTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method = JetTagging::JetTaggingSFMethod::mujets, const variation syst = variation::nom, const TString &source = "total");
    float GetBTaggingR(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, std::string &processName, const variation syst = variation::nom, const TString &source = "total") const;
    inline float GetBTaggingSF(const RVec<Jet> &jets, const JetTagging::JetTaggingSFMethod method = JetTagging::JetTaggingSFMethod::mujets, const variation syst = variation::nom, const TString &source = "total") { return GetBTaggingSF(jets, global_tagger, global_wp, method, syst, source); }
    inline float GetBTaggingR(const RVec<Jet> &jets, std::string &processName, const variation syst = variation::nom, const TString &source = "total") const { return GetBTaggingR(jets, global_tagger, processName, syst, source); }

    // ctagging
    pair<float, float> GetCTaggingWP() const;
    pair<float, float> GetCTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const;
    float GetCTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst = variation::nom);
    float GetCTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method, const variation syst, const TString &source = "total");
    float GetCTaggingR(const float npvs, const float HT, const JetTagging::JetFlavTagger tagger, const TString &processName = "", const TString &ttBarCategory = "total", const TString &syst_str = "") const;
    inline float GetCTaggingSF(const RVec<Jet> &jets, const JetTagging::JetTaggingSFMethod method = JetTagging::JetTaggingSFMethod::mujets, const variation syst = variation::nom, const TString &source = "total") { return GetCTaggingSF(jets, global_tagger, global_wp, method, syst, source); }

    // PileUp Jet ID
    float GetPileupJetIDSF(const RVec<Jet> &jets, const unordered_map<int, int> &matched_idx, const TString &wp, const variation syst=variation::nom);

    // JERC
    float GetJER(const float eta, const float pt, const float rho) const;
    float GetJERSF(const float eta, const float pt, const variation syst = variation::nom, const TString &source = "total") const;
    float GetJESSF(const float area, const float eta, const float pt, const float phi, const float rho, const unsigned int runNumber) const;
    float GetJESUncertainty(const float eta, const float pt, const variation syst = variation::nom, const TString &source = "total") const;
    // jerc_fatjet
    
    // jetvetomap
    bool IsJetVetoZone(const float eta, const float phi, TString mapCategory) const;
    
    // MET
    void METXYCorrection(Particle &Met, const int RunNumber, const int npvs, const XYCorrection_MetType MetType);
    
    
    // reweighting
    float GetTopPtReweight(const RVec<Gen> &gens) const;

    // WZ Nj SF
    float GetWZNjetsSF(const float &njets, const string &syst="Central") const;

    // Safe evaluation function for correction sets with comprehensive error handling
    template<typename... Args>
    inline float safeEvaluate(const correction::Correction::Ref &cset, 
                              const string &function_name,
                              const vector<correction::Variable::Type> &args) const {
        if (!cset) {
            cerr << "[MyCorrection::" << function_name << "] Error: Correction set is null" << endl;
            exit(EXIT_FAILURE);
        }
        
        try {
            return cset->evaluate(args);
        } catch (const std::exception &e) {
            cerr << "[MyCorrection::" << function_name << "] Error during evaluation: " << e.what() << endl;
            cerr << "[MyCorrection::" << function_name << "] Arguments (" << args.size() << "): ";
            for (const auto &arg : args) {
                std::visit([](const auto &value) { cerr << value << " "; }, arg);
            }
            cerr << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Overload for CompoundCorrection
    template<typename... Args>
    inline float safeEvaluate(const correction::CompoundCorrection::Ref &cset, 
                              const string &function_name,
                              const vector<correction::Variable::Type> &args) const {
        if (!cset) {
            cerr << "[MyCorrection::" << function_name << "] Error: CompoundCorrection set is null" << endl;
            exit(EXIT_FAILURE);
        }
        
        try {
            return cset->evaluate(args);
        } catch (const std::exception &e) {
            cerr << "[MyCorrection::" << function_name << "] Error during evaluation: " << e.what() << endl;
            cerr << "[MyCorrection::" << function_name << "] Arguments (" << args.size() << "): ";
            for (const auto &arg : args) {
                std::visit([](const auto &value) { cerr << value << " "; }, arg);
            }
            cerr << endl;
            exit(EXIT_FAILURE);
        }
    }

private:
    struct EraConfig {
        string json_muon;
        string json_muon_trig_eff;
        string json_puWeights;
        string json_btagging;
        string json_ctagging;
        string json_btagging_eff;
        string json_ctagging_eff;
        string json_btagging_R;
        string json_ctagging_R;
        string json_electron;
        string json_electron_hlt;
        string json_electron_variation;
        string json_photon;
        string json_jerc;
        string json_jerc_fatjet;
        string json_jetvetomap;
        string json_jmar;
        string json_met;
        string txt_roccor;
        
        // Custom
        string json_muon_TopHNT_idsf;
        string json_muon_TopHNT_dblmu_leg1_eff;
        string json_muon_TopHNT_dblmu_leg2_eff;
        string json_muon_TopHNT_emu_leg1_eff;
        string json_muon_TopHNT_emu_leg2_eff;
        string json_muon_TopHNT_fakerate;
        string json_electron_TopHNT_idsf;
        string json_electron_TopHNT_emu_leg1_eff;
        string json_electron_TopHNT_emu_leg2_eff;
        string json_electron_TopHNT_fakerate;
        string json_WZ_Nj_SF;
    };
    EraConfig GetEraConfig(TString era, const string &btagging_eff_file, const string &ctagging_eff_file, const string &btagging_R_file, const string &ctagging_R_file) const;

    inline void SetEra(TString era) {
        DataEra = era;
        if (era.Contains("2016") || era.Contains("2017") || era.Contains("2018")) {
            Run = 2;
        } else if (era.Contains("2022") || era.Contains("2023")) {
            Run = 3;
        } else {
            throw runtime_error("Invalid era: " + era);
        }
    }
    inline void SetPeriod(TString period) { DataPeriod = period; }
    inline TString GetEra() const { return DataEra; }
    inline TString GetPeriod() const { return DataPeriod; }
    inline void SetSample(TString sample) { Sample = sample; }
    inline void setIsData(bool isData) { IsDATA = isData; }

    inline bool loadCorrectionSet(const string &name,
                                  const string &file,
                                  unique_ptr<CorrectionSet> &cset,
                                  bool optional = false) {
        cout << "[MyCorrection::loadCorrectionSet] " << name << ": " << file << endl;
        try {
            cset = CorrectionSet::from_file(file);
            return true;
        } catch (const exception &e) {
            if (optional) {
                cerr << "[MyCorrection::loadCorrectionSet] Warning: Failed to load " << name << " (" << file << "): " << e.what() << endl;
                return false;
            } else {
                cerr << "[MyCorrection::loadCorrectionSet] Error: Failed to load " << name << " (" << file << "): " << e.what() << endl;
                throw;
            }
        }
    }

    inline bool loadRoccoR(const string &file, bool optional=false) {
        cout << "[MyCorrection::loadRoccoR] " << file << endl;
        try {
            rc.init(file);
            return true;
        } catch (const exception &e) {
            if (optional) {
                cerr << "[MyCorrection::loadRoccoR] Warning: Failed to load " << file << " (" << e.what() << ")" << endl;
                return false;
            } else {
                cerr << "[MyCorrection::loadRoccoR] Error: Failed to load " << file << " (" << e.what() << ")" << endl;
                throw;
            }
        }
    }

    inline bool isInputInCorrection(const string &key, const correction::Correction::Ref &cset) const {
        vector<string> inputs;
        for (const auto &input : cset->inputs()) {
            inputs.push_back(input.name());
        }
        return find(inputs.begin(), inputs.end(), key) != inputs.end();
    }

    JetTagging::JetFlavTaggerWP global_wp;
    JetTagging::JetFlavTagger global_tagger;
    string global_wpStr;
    string global_taggerStr;
    TString DataEra;
    TString DataPeriod;
    int Run;
    TString Sample;
    bool IsDATA;

    unique_ptr<CorrectionSet> cset_muon;
    unique_ptr<CorrectionSet> cset_muon_trig_eff;
    unique_ptr<CorrectionSet> cset_puWeights;
    unique_ptr<CorrectionSet> cset_btagging;
    unique_ptr<CorrectionSet> cset_ctagging;
    unique_ptr<CorrectionSet> cset_btagging_eff;
    unique_ptr<CorrectionSet> cset_ctagging_eff;
    unique_ptr<CorrectionSet> cset_btagging_R;
    unique_ptr<CorrectionSet> cset_ctagging_R;
    unique_ptr<CorrectionSet> cset_electron;
    unique_ptr<CorrectionSet> cset_electron_hlt;
    unique_ptr<CorrectionSet> cset_electron_variation;
    unique_ptr<CorrectionSet> cset_photon;
    unique_ptr<CorrectionSet> cset_jerc;
    unique_ptr<CorrectionSet> cset_jerc_fatjet;
    unique_ptr<CorrectionSet> cset_jetvetomap;
    unique_ptr<CorrectionSet> cset_jmar;
    unique_ptr<CorrectionSet> cset_met;

    // custom
    unique_ptr<CorrectionSet> cset_muon_TopHNT_idsf;
    unique_ptr<CorrectionSet> cset_muon_TopHNT_dblmu_leg1_eff;
    unique_ptr<CorrectionSet> cset_muon_TopHNT_dblmu_leg2_eff;
    unique_ptr<CorrectionSet> cset_muon_TopHNT_emu_leg1_eff;
    unique_ptr<CorrectionSet> cset_muon_TopHNT_emu_leg2_eff;
    unique_ptr<CorrectionSet> cset_muon_TopHNT_fakerate;
    unique_ptr<CorrectionSet> cset_electron_TopHNT_idsf;
    unique_ptr<CorrectionSet> cset_electron_TopHNT_emu_leg1_eff;
    unique_ptr<CorrectionSet> cset_electron_TopHNT_emu_leg2_eff;
    unique_ptr<CorrectionSet> cset_electron_TopHNT_fakerate;
    unique_ptr<CorrectionSet> cset_WZ_Nj_SF;

    unordered_map<string, string> MUO_keys;
    unordered_map<string, string> LUM_keys;
    unordered_map<string, string> BTV_keys;
    unordered_map<string, string> EGM_keys;
    unordered_map<string, string> JME_JER_GT;
    unordered_map<string, string> JME_JES_GT;
    unordered_map<string, string> JME_vetomap_keys;
    unordered_map<string, string> JME_PILEUP_keys;
    unordered_map<string, string> JME_MET_keys;
    
    RoccoR rc;

    // All POG choose different string for the systematics, so we need to convert the enum to the string....
    // Here I implement every single function instead of a general one, because heavy use of switch-case might be slow.
    inline string getSystString_CUSTOM(const variation syst) const {
        string sys_string = "nom";
        switch (syst) {
        case variation::nom:
            sys_string = "nom";
            break;
        case variation::up:
            sys_string = "up";
            break;
        case variation::down:
            sys_string = "down";
            break;
        };
        return sys_string;
    };

    inline string getSystString_MUO(const variation syst) const {
        string sys_string = "nominal";
        switch (syst) {
        case variation::nom:
            sys_string = "nominal";
            break;
        case variation::up:
            sys_string = "systup";
            break;
        case variation::down:
            sys_string = "systdown";
            break;
        };
        return sys_string;
    };

    inline string getSystString_LUM(const variation syst) const {
        string sys_string = "nominal";
        switch (syst) {
        case variation::nom:
            sys_string = "nominal";
            break;
        case variation::up:
            sys_string = "up";
            break;
        case variation::down:
            sys_string = "down";
            break;
        };
        return sys_string;
    };

    inline string getSystString_BTV(const variation syst) const {
        string sys_string = "central";
        switch (syst) {
        case variation::nom:
            sys_string = "central";
            break;
        case variation::up:
            sys_string = "up";
            break;
        case variation::down:
            sys_string = "down";
            break;
        };
        return sys_string;
    };

    inline string getSystString_EGM(const variation syst) const {
        string sys_string = "sf";
        switch (syst) {
        case variation::nom:
            sys_string = "sf";
            break;
        case variation::up:
            sys_string = "sfup";
            break;
        case variation::down:
            sys_string = "sfdown";
            break;
        };
        return sys_string;
    };
    
    inline string getSystString_EGMScale(const variation syst) const {
        // Only for Run2
        if (! (Run == 2)) {
            throw runtime_error("[MyCorrection::getSystString_EGMScale] Use getSystString_EGM for Run3");
        }
        string sys_string = "";
        switch(syst) {
        case variation::nom:
            sys_string = "";
            break;
        case variation::up:
            sys_string = "scaleup";
            break;
        case variation::down:
            sys_string = "scaledown";
            break;
        };
        return sys_string;
    }
    
    inline string getSystString_JME(const variation syst) const {
        string sys_string = "nom";
        switch (syst) {
        case variation::nom:
            sys_string = "nom";
            break;
        case variation::up:
            sys_string = "up";
            break;
        case variation::down:
            sys_string = "down";
            break;
        };
        return sys_string;
    };
};

#endif
