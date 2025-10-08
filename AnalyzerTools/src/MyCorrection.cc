#include "MyCorrection.h"

MyCorrection::MyCorrection() {}
MyCorrection::MyCorrection(const TString &era, const TString &period, const TString &sample, const bool IsData, const string &btagging_eff_file, const string &ctagging_eff_file, const string &btagging_R_file, const string &ctagging_R_file) { 
    cout << "[MyCorrection::MyCorrection] MyCorrection created for " << era << endl;
    SetEra(era);
    SetPeriod(period);
    SetSample(sample);
    setIsData(IsData);

    EraConfig config = GetEraConfig(era, btagging_eff_file, ctagging_eff_file, btagging_R_file, ctagging_R_file);
    struct CorrectionInfo {
        string name;
        string path;  
        unique_ptr<CorrectionSet> &cset;
        bool isOptional;
    };

    // load correction sets
    loadCorrectionSet("muon SF", config.json_muon, cset_muon, false);
    loadCorrectionSet("puWeights", config.json_puWeights, cset_puWeights, false);
    loadCorrectionSet("btagging", config.json_btagging, cset_btagging, false);
    loadCorrectionSet("ctagging", config.json_ctagging, cset_ctagging, false);
    loadCorrectionSet("btagging eff", config.json_btagging_eff, cset_btagging_eff, false);
    loadCorrectionSet("ctagging eff", config.json_ctagging_eff, cset_ctagging_eff, false);
    loadCorrectionSet("electron", config.json_electron, cset_electron, false);
    loadCorrectionSet("electron variation", config.json_electron_variation, cset_electron_variation, false);
    loadCorrectionSet("photon", config.json_photon, cset_photon, false);
    loadCorrectionSet("jerc", config.json_jerc, cset_jerc, false);
    loadCorrectionSet("jerc_fatjet", config.json_jerc_fatjet, cset_jerc_fatjet, false);
    loadCorrectionSet("jetvetomap", config.json_jetvetomap, cset_jetvetomap, false);
    // Optional files
    loadRoccoR(config.txt_roccor, true);
    loadCorrectionSet("jmar", config.json_jmar, cset_jmar, true);
    loadCorrectionSet("muon trig eff", config.json_muon_trig_eff, cset_muon_trig_eff, true);
    loadCorrectionSet("electron hlt", config.json_electron_hlt, cset_electron_hlt, true);
    loadCorrectionSet("met", config.json_met, cset_met, true);
    loadCorrectionSet("btagging R", config.json_btagging_R, cset_btagging_R, true);
    loadCorrectionSet("ctagging R", config.json_ctagging_R, cset_ctagging_R, true);
    loadCorrectionSet("muon TopHNT idsf", config.json_muon_TopHNT_idsf, cset_muon_TopHNT_idsf, true);
    loadCorrectionSet("muon TopHNT dblmu leg1 eff", config.json_muon_TopHNT_dblmu_leg1_eff, cset_muon_TopHNT_dblmu_leg1_eff, true);
    loadCorrectionSet("muon TopHNT dblmu leg2 eff", config.json_muon_TopHNT_dblmu_leg2_eff, cset_muon_TopHNT_dblmu_leg2_eff, true);
    loadCorrectionSet("muon TopHNT emu leg1 eff", config.json_muon_TopHNT_emu_leg1_eff, cset_muon_TopHNT_emu_leg1_eff, true);
    loadCorrectionSet("muon TopHNT emu leg2 eff", config.json_muon_TopHNT_emu_leg2_eff, cset_muon_TopHNT_emu_leg2_eff, true);
    loadCorrectionSet("muon TopHNT fakerate", config.json_muon_TopHNT_fakerate, cset_muon_TopHNT_fakerate, true);
    loadCorrectionSet("electron TopHNT idsf", config.json_electron_TopHNT_idsf, cset_electron_TopHNT_idsf, true);
    loadCorrectionSet("electron TopHNT emu leg1 eff", config.json_electron_TopHNT_emu_leg1_eff, cset_electron_TopHNT_emu_leg1_eff, true);
    loadCorrectionSet("electron TopHNT emu leg2 eff", config.json_electron_TopHNT_emu_leg2_eff, cset_electron_TopHNT_emu_leg2_eff, true);
    loadCorrectionSet("electron TopHNT fakerate", config.json_electron_TopHNT_fakerate, cset_electron_TopHNT_fakerate, true);
    loadCorrectionSet("WZ Nj SF", config.json_WZ_Nj_SF, cset_WZ_Nj_SF, true);


    LUM_keys["2023BPix"] = "Collisions2023_369803_370790_eraD_GoldenJson";
    LUM_keys["2023"] = "Collisions2023_366403_369802_eraBC_GoldenJson";
    LUM_keys["2022EE"] = "Collisions2022_359022_362760_eraEFG_GoldenJson";
    LUM_keys["2022"] = "Collisions2022_355100_357900_eraBCD_GoldenJson";
    LUM_keys["2018"] = "Collisions18_UltraLegacy_goldenJSON";
    LUM_keys["2017"] = "Collisions17_UltraLegacy_goldenJSON";
    LUM_keys["2016postVFP"] = "Collisions16_UltraLegacy_goldenJSON";
    LUM_keys["2016preVFP"] = "Collisions16_UltraLegacy_goldenJSON";

    EGM_keys["2023BPix"] = "2023PromptD";
    EGM_keys["2023"] = "2023PromptC";
    EGM_keys["2022EE"] = "2022Re-recoE+PromptFG";
    EGM_keys["2022"] = "2022Re-recoBCD";
    EGM_keys["2016preVFP"] = "2016preVFP";
    EGM_keys["2016postVFP"] = "2016postVFP";
    EGM_keys["2017"] = "2017";
    EGM_keys["2018"] = "2018";

    // JER tags
    // ###### as placeholder
    JME_JER_GT["2023BPix"] = "Summer23BPixPrompt23_RunD_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2023"] = "Summer23Prompt23_RunCv1234_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2022EE"] = "Summer22EE_22Sep2023_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2022"] = "Summer22_22Sep2023_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2018"] = "Summer19UL18_JRV2_MC_######_AK4PFchs";
    JME_JER_GT["2017"] = "Summer19UL17_JRV2_MC_######_AK4PFchs";
    JME_JER_GT["2016postVFP"] = "Summer20UL16_JRV3_MC_######_AK4PFchs";
    JME_JER_GT["2016preVFP"] = "Summer20UL16APV_JRV3_MC_######_AK4PFchs";
    
    // JES tags
    // ###### as placeholder
    if (IsDATA) {
        string period = static_cast<string>(GetPeriod().Data());
        if (GetEra() == "2016preVFP") {
            period = TString("EF").Contains(GetPeriod()) ? "EF" : "BCD";
        }
        JME_JES_GT["2023BPix"] = "Summer23BPixPrompt23_V3_DATA_######_AK4PFPuppi";
        JME_JES_GT["2023"] = "Summer23Prompt23_V2_DATA_######_AK4PFPuppi";
        JME_JES_GT["2022EE"] = "Summer22EE_22Sep2023_Run"+period+"_V2_DATA_######_AK4PFPuppi";
        JME_JES_GT["2022"] = "Summer22_22Sep2023_RunCD_V2_DATA_######_AK4PFPuppi";
        JME_JES_GT["2018"] = "Summer19UL18_Run"+period+"_V5_DATA_######_AK4PFchs";
        JME_JES_GT["2017"] = "Summer19UL17_Run"+period+"_V5_DATA_######_AK4PFchs";
        JME_JES_GT["2016postVFP"] = "Summer19UL16_RunFGH_V7_DATA_######_AK4PFchs";
        JME_JES_GT["2016preVFP"] = "Summer19UL16APV_Run"+period+"_V7_DATA_######_AK4PFchs";
    } else { // MC
        JME_JES_GT["2023BPix"] = "Summer23BPixPrompt23_V3_MC_######_AK4PFPuppi";
        JME_JES_GT["2023"] = "Summer23Prompt23_V2_MC_######_AK4PFPuppi";
        JME_JES_GT["2022EE"] = "Summer22EE_22Sep2023_V2_MC_######_AK4PFPuppi";
        JME_JES_GT["2022"] = "Summer22_22Sep2023_V2_MC_######_AK4PFPuppi";
        JME_JES_GT["2018"] = "Summer19UL18_V5_MC_######_AK4PFchs";
        JME_JES_GT["2017"] = "Summer19UL17_V5_MC_######_AK4PFchs";
        JME_JES_GT["2016postVFP"] = "Summer19UL16_V7_MC_######_AK4PFchs";
        JME_JES_GT["2016preVFP"] = "Summer19UL16APV_V7_MC_######_AK4PFchs"; 
    }
    
    // Jet VetpMap
    JME_vetomap_keys["2023BPix"] = "Summer23BPixPrompt23_RunD_V1";
    JME_vetomap_keys["2023"] = "Summer23Prompt23_RunC_V1";
    JME_vetomap_keys["2022EE"] = "Summer22EE_23Sep2023_RunEFG_V1";
    JME_vetomap_keys["2022"] = "Summer22_23Sep2023_RunCD_V1";
    JME_vetomap_keys["2018"] = "Summer19UL18_V1";
    JME_vetomap_keys["2017"] = "Summer19UL17_V1";
    JME_vetomap_keys["2016postVFP"] = "Summer19UL16_V1";
    JME_vetomap_keys["2016preVFP"] = "Summer19UL16_V1";

    JME_PILEUP_keys["2016preVFP"] = "PUJetID_eff";
    JME_PILEUP_keys["2016postVFP"] = "PUJetID_eff";
    JME_PILEUP_keys["2017"] = "PUJetID_eff";
    JME_PILEUP_keys["2018"] = "PUJetID_eff";
}

MyCorrection::~MyCorrection() {}

MyCorrection::EraConfig MyCorrection::GetEraConfig(TString era, const string &btagging_eff_file, const string &ctagging_eff_file, const string &btagging_R_file, const string &ctagging_R_file) const {
    EraConfig config;

    const char *json_pog_path = getenv("JSONPOG_REPO_PATH");
    const char *sknano_data = getenv("SKNANO_DATA");
    const char *external_roccor = getenv("ROCCOR_PATH");
    cout << "[MyCorrection::GetEraConfig] json_pog_path: " << (json_pog_path ? json_pog_path : "NULL") << endl;
    cout << "[MyCorrection::GetEraConfig] sknano_data: " << (sknano_data ? sknano_data : "NULL") << endl;
    cout << "[MyCorrection::GetEraConfig] external_roccor: " << (external_roccor ? external_roccor : "NULL") << endl;
    if (!json_pog_path || !sknano_data || !external_roccor) {
        throw runtime_error("JSONPOG_REPO_PATH or SKNANO_DATA or ROCCOR_PATH is not set");
    }

    const string json_pog_path_str(json_pog_path);
    const string sknano_data_str(sknano_data);
    const string external_roccor_str(external_roccor);

    config.json_muon = json_pog_path_str + "/POG/MUO";
    config.json_muon_trig_eff = sknano_data_str;
    config.json_puWeights = json_pog_path_str + "/POG/LUM";
    config.json_btagging = json_pog_path_str + "/POG/BTV";
    config.json_ctagging = json_pog_path_str + "/POG/BTV";
    config.json_btagging_eff = sknano_data_str;
    config.json_ctagging_eff = sknano_data_str;
    config.json_btagging_R = sknano_data_str;
    config.json_ctagging_R = sknano_data_str;
    config.json_electron = json_pog_path_str + "/POG/EGM";
    config.json_electron_hlt = config.json_electron;
    config.json_electron_variation = sknano_data_str;
    config.json_photon = json_pog_path_str + "/POG/EGM";
    config.json_jerc = json_pog_path_str + "/POG/JME";
    config.json_jerc_fatjet = json_pog_path_str + "/POG/JME";
    config.json_jetvetomap = json_pog_path_str + "/POG/JME";
    config.json_jmar = json_pog_path_str + "/POG/JME";
    config.json_met = json_pog_path_str + "/POG/JME";
    config.txt_roccor = external_roccor_str;

    config.json_muon_TopHNT_idsf = sknano_data_str + "/" + DataEra.Data() + "/MUO/efficiency_TopHNT.json";
    config.json_muon_TopHNT_dblmu_leg1_eff = sknano_data_str + "/" + DataEra.Data() + "/MUO/efficiency_Mu17Leg1.json";
    config.json_muon_TopHNT_dblmu_leg2_eff = sknano_data_str + "/" + DataEra.Data() + "/MUO/efficiency_Mu8Leg2.json";
    config.json_muon_TopHNT_emu_leg1_eff = sknano_data_str + "/" + DataEra.Data() + "/MUO/efficiency_Mu23El12_Mu23Leg.json";
    config.json_muon_TopHNT_emu_leg2_eff = sknano_data_str + "/" + DataEra.Data() + "/MUO/efficiency_Mu8El23_Mu8Leg.json";
    config.json_muon_TopHNT_fakerate = sknano_data_str + "/" + DataEra.Data() + "/MUO/fakerate_TopHNT.json";
    config.json_electron_TopHNT_idsf = sknano_data_str + "/" + DataEra.Data() + "/EGM/efficiency_TopHNT.json";
    config.json_electron_TopHNT_emu_leg1_eff = sknano_data_str + "/" + DataEra.Data() + "/EGM/efficiency_Mu8El23_El23Leg.json";
    config.json_electron_TopHNT_emu_leg2_eff = sknano_data_str + "/" + DataEra.Data() + "/EGM/efficiency_Mu23El12_El12Leg.json";
    config.json_electron_TopHNT_fakerate = sknano_data_str + "/" + DataEra.Data() + "/EGM/fakerate_TopHNT.json";
    config.json_WZ_Nj_SF = sknano_data_str;

    if (era == "2016preVFP") {
        config.json_muon += "/2016preVFP_UL/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2016preVFP/MUO/muon_trig.json";
        config.json_puWeights += "/2016preVFP_UL/puWeights.json.gz";
        config.json_btagging += "/2016preVFP_UL/btagging.json.gz";
        config.json_ctagging += "/2016preVFP_UL/ctagging.json.gz";
        config.json_btagging_eff += "/2016preVFP/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2016preVFP/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2016preVFP/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2016preVFP/BTV/" + ctagging_R_file;
        config.json_electron += "/2016preVFP_UL/electron.json.gz";
        config.json_electron_variation += "/2016preVFP/EGM/EGM_ScaleUnc.json.gz";
        config.json_photon += "/2016preVFP_UL/photon.json.gz";
        config.json_jerc += "/2016preVFP_UL/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2016preVFP_UL/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2016preVFP_UL/jetvetomaps.json.gz";
        config.json_jmar += "/2016preVFP_UL/jmar.json.gz";
        config.json_met += "/2016preVFP_UL/met.json.gz";
        config.txt_roccor += "/RoccoR2016aUL.txt";
    } else if (era == "2016postVFP") {
        config.json_muon += "/2016postVFP_UL/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2016postVFP/MUO/muon_trig.json";
        config.json_puWeights += "/2016postVFP_UL/puWeights.json.gz";
        config.json_btagging += "/2016postVFP_UL/btagging.json.gz";
        config.json_ctagging += "/2016postVFP_UL/ctagging.json.gz";
        config.json_btagging_eff += "/2016postVFP/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2016postVFP/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2016postVFP/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2016postVFP/BTV/" + ctagging_R_file;
        config.json_electron += "/2016postVFP_UL/electron.json.gz";
        config.json_electron_variation += "/2016postVFP/EGM/EGM_ScaleUnc.json.gz";
        config.json_photon += "/2016postVFP_UL/photon.json.gz";
        config.json_jerc += "/2016postVFP_UL/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2016postVFP_UL/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2016postVFP_UL/jetvetomaps.json.gz";
        config.json_jmar += "/2016postVFP_UL/jmar.json.gz";
        config.json_met += "/2016postVFP_UL/met.json.gz";
        config.txt_roccor += "/RoccoR2016bUL.txt";
    } else if (era == "2017") {
        config.json_muon += "/2017_UL/muon_Z.json.gz";
        //config.json_muon_trig_eff += "/2017/MUO/muon_trig.json";
        config.json_puWeights += "/2017_UL/puWeights.json.gz";
        config.json_btagging += "/2017_UL/btagging.json.gz";
        config.json_ctagging += "/2017_UL/ctagging.json.gz";
        config.json_btagging_eff += "/2017/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2017/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2017/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2017/BTV/" + ctagging_R_file;
        config.json_electron += "/2017_UL/electron.json.gz";
        config.json_electron_variation += "/2017/EGM/EGM_ScaleUnc.json.gz";
        config.json_photon += "/2017_UL/photon.json.gz";
        config.json_jerc += "/2017_UL/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2017_UL/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2017_UL/jetvetomaps.json.gz";
        config.json_jmar += "/2017_UL/jmar.json.gz";
        config.json_met += "/2017_UL/met.json.gz";
        config.txt_roccor += "/RoccoR2017UL.txt";
    } else if (era == "2018") {
        config.json_muon += "/2018_UL/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2018/MUO/muon_trig.json";
        config.json_puWeights += "/2018_UL/puWeights.json.gz";
        config.json_btagging += "/2018_UL/btagging.json.gz";
        config.json_ctagging += "/2018_UL/ctagging.json.gz";
        config.json_btagging_eff += "/2018/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2018/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2018/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2018/BTV/" + ctagging_R_file;
        config.json_electron += "/2018_UL/electron.json.gz";
        config.json_electron_variation += "/2018/EGM/EGM_ScaleUnc.json.gz";
        config.json_photon += "/2018_UL/photon.json.gz";
        config.json_jerc += "/2018_UL/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2018_UL/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2018_UL/jetvetomaps.json.gz";
        config.json_jmar += "/2018_UL/jmar.json.gz";
        config.json_met += "/2018_UL/met.json.gz";
        config.txt_roccor += "/RoccoR2018UL.txt";
    } else if (era == "2022") {
        config.json_muon += "/2022_Summer22/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2022/MUO/muon_trig.json";
        config.json_puWeights += "/2022_Summer22/puWeights.json.gz";
        config.json_btagging += "/2022_Summer22/btagging.json.gz";
        config.json_ctagging += "/2022_Summer22/ctagging.json.gz";
        config.json_btagging_eff += "/2022/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2022/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2022/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2022/BTV/" + ctagging_R_file;
        config.json_electron += "/2022_Summer22/electron.json.gz";
        config.json_electron_variation = json_pog_path_str + "/POG/EGM/2022_Summer22/electronSS.json.gz";
        config.json_electron_hlt += "/2022_Summer22/electronHlt.json.gz";
        config.json_photon += "/2022_Summer22/photon.json.gz";
        config.json_jerc += "/2022_Summer22/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2022_Summer22/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2022_Summer22/jetvetomaps.json.gz";
        config.json_met += "/2022_Summer22/met.json.gz";
        config.txt_roccor += "/RoccoR2022.txt";
    } else if (era == "2022EE") {
        config.json_muon += "/2022_Summer22EE/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2022EE/MUO/muon_trig.json";
        config.json_puWeights += "/2022_Summer22EE/puWeights.json.gz";
        config.json_btagging += "/2022_Summer22EE/btagging.json.gz";
        config.json_ctagging += "/2022_Summer22EE/ctagging.json.gz";
        config.json_btagging_eff += "/2022EE/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2022EE/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2022EE/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2022EE/BTV/" + ctagging_R_file;
        config.json_electron += "/2022_Summer22EE/electron.json.gz";
        config.json_electron_variation = json_pog_path_str + "/POG/EGM/2022_Summer22EE/electronSS.json.gz";
        config.json_electron_hlt += "/2022_Summer22EE/electronHlt.json.gz";
        config.json_photon += "/2022_Summer22EE/photon.json.gz";
        config.json_jerc += "/2022_Summer22EE/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2022_Summer22EE/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2022_Summer22EE/jetvetomaps.json.gz";
        config.json_met += "/2022_Summer22EE/met.json.gz";
        config.txt_roccor += "/RoccoR2022EE.txt";
    } else if (era == "2023") {
        config.json_muon += "/2023_Summer23/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2023/MUO/muon_trig.json";
        config.json_puWeights += "/2023_Summer23/puWeights.json.gz";
        config.json_btagging += "/2023_Summer23/btagging.json.gz";
        config.json_ctagging += "/2023_Summer23/ctagging.json.gz";
        config.json_btagging_eff += "/2023/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2023/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2023/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2023/BTV/" + ctagging_R_file;
        config.json_electron += "/2023_Summer23/electron.json.gz";
        config.json_electron_variation = json_pog_path_str + "/POG/EGM/2023_Summer23/electronSS.json.gz";
        config.json_electron_hlt += "/2023_Summer23/electronHlt.json.gz";
        config.json_photon += "/2023_Summer23/photon.json.gz";
        config.json_jerc += "/2023_Summer23/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2023_Summer23/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2023_Summer23/jetvetomaps.json.gz";
        config.json_met += "/2023_Summer23/met.json.gz";
        config.txt_roccor += "/RoccoR2023.txt";
    } else if (era == "2023BPix") {
        config.json_muon += "/2023_Summer23BPix/muon_Z.json.gz";
        config.json_muon_trig_eff += "/2023BPix/MUO/muon_trig.json";
        config.json_puWeights += "/2023_Summer23BPix/puWeights.json.gz";
        config.json_btagging += "/2023_Summer23BPix/btagging.json.gz";
        config.json_ctagging += "/2023_Summer23BPix/ctagging.json.gz";
        config.json_btagging_eff += "/2023BPix/BTV/" + btagging_eff_file;
        config.json_ctagging_eff += "/2023BPix/BTV/" + ctagging_eff_file;
        config.json_btagging_R += "/2023BPix/BTV/" + btagging_R_file;
        config.json_ctagging_R += "/2023BPix/BTV/" + ctagging_R_file;
        config.json_electron += "/2023_Summer23BPix/electron.json.gz";
        config.json_electron_variation = json_pog_path_str + "/POG/EGM/2023_Summer23BPix/electronSS.json.gz";
        config.json_electron_hlt += "/2023_Summer23BPix/electronHlt.json.gz";
        config.json_photon += "/2023_Summer23BPix/photon.json.gz";
        config.json_jerc += "/2023_Summer23BPix/jet_jerc.json.gz";
        config.json_jerc_fatjet += "/2023_Summer23BPix/fatJet_jerc.json.gz";
        config.json_jetvetomap += "/2023_Summer23BPix/jetvetomaps.json.gz";
        config.json_met += "/2023_Summer23BPix/met.json.gz";
        config.txt_roccor += "/RoccoR2023BPix.txt";
    } else {
        throw invalid_argument("[MyCorrection::GetEraConfig] Invalid era: " + era);
    }

    if (Run == 2) {
        config.json_WZ_Nj_SF += "/Run2/WZSF/WZNjetsSF.json";
    } else if (Run == 3) {
        config.json_WZ_Nj_SF += "/Run3/WZSF/WZNjetsSF.json";
    } else {
        throw invalid_argument("[MyCorrection::GetEraConfig] Invalid run number");
    }
    
    return config;
}

// Muon
// Rochestor correction
float MyCorrection::GetMuonScaleSF(const Muon &muon, const variation syst, const float matched_pt) const {
    float roccor = 1.;
    float roccor_err = 0.;

    // few GeVs of muon shuold not use this correction, because the authors did not consider the radiations of low pt muons inside detectors
    // still true for Run3?
    if (muon.Pt() < 10.) return 1.;

    if (IsDATA) {
        roccor = rc.kScaleDT(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), 0, 0);
        roccor_err = rc.kScaleDTerror(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi());
    } else {
        // Random seed is initialized in SKNanoLoader::Init()
        gRandom->SetSeed(int(muon.Pt()/muon.Eta()));
        float u = gRandom->Rndm();
        if (matched_pt > 0) { // matched
            roccor = rc.kSpreadMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), matched_pt, 0, 0);
            roccor_err = rc.kSpreadMCerror(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), matched_pt);
        } else {
            //roccor = rc.kScaleMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), 0, 0);
            //roccor_err = 0.;
            //roccor_err = rc.kScaleMCerror(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi());
            roccor = rc.kSmearMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), muon.nTrackerLayers(), u, 0, 0);
            roccor_err = rc.kSmearMCerror(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), muon.nTrackerLayers(), u);
        }
    }

    if (syst == variation::nom) return roccor;
    else if (syst == variation::up) return roccor + roccor_err;
    else if (syst == variation::down) return roccor - roccor_err;
    else {
        throw runtime_error("[MyCorrection::GetMuonScaleSF] Invalid syst value");
    }
}

float MyCorrection::GetMuonRECOSF(const Muon &muon, const variation syst) const {
    // No correction for Run3, see https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_5_2
    if (Run == 3) return 1.;

    // For RECO efficiency, used 40-60 GeV muons due to the large background in Z-peak.
    // Plaetue already reached in a few GeV, okay to use for [10, 200] GeV muons. 
    auto cset = cset_muon->at("NUM_TrackerMuons_DEN_genTracks");
    return safeEvaluate(cset, "GetMuonRECOSF", {muon.Eta(), (muon.OriginalPt() < 40. ? 40. : muon.OriginalPt()), getSystString_MUO(syst)});
}

float MyCorrection::GetMuonRECOSF(const RVec<Muon> &muons, const variation syst) const {
    float weight = 1.;
    for (const auto &muon : muons) {
        weight *= GetMuonRECOSF(muon, syst);
    }
    return weight;
}

float MyCorrection::GetMuonIDSF(const TString &Muon_ID_SF_Key, const Muon &muon, const variation syst) const {
    if (Muon_ID_SF_Key == "TopHNT") {
        auto cset = cset_muon_TopHNT_idsf->at("sf");
        if (syst == variation::nom) {
            return safeEvaluate(cset, "GetMuonIDSF", {fabs(muon.Eta()), muon.OriginalPt(), "nom"});
        } else if (syst == variation::up) {
            return safeEvaluate(cset, "GetMuonIDSF", {fabs(muon.Eta()), muon.OriginalPt(), "up"});
        } else if (syst == variation::down) {
            return safeEvaluate(cset, "GetMuonIDSF", {fabs(muon.Eta()), muon.OriginalPt(), "down"});
        } else {
            throw runtime_error("[MyCorrection::GetMuonIDSF] Invalid syst value");
        }
    } else {
        auto cset = cset_muon->at(string(Muon_ID_SF_Key));
        return safeEvaluate(cset, "GetMuonIDSF", {fabs(muon.Eta()), muon.OriginalPt(), getSystString_MUO(syst)});
    }
}

float MyCorrection::GetMuonIDSF(const TString &Muon_ID_SF_Key, const RVec<Muon> &muons, const variation syst) const {
    float weight = 1.;
    for (const auto &muon : muons) {
        weight *= GetMuonIDSF(Muon_ID_SF_Key, muon, syst);
    }
    return weight;
}

float MyCorrection::GetFakeRate(const Muon muon, const TString &id_key, const TString &syst_key) const {
    if (id_key == "TopHNT") {
        const float absEta = fabs(muon.Eta());
        float ptCorr = muon.Pt()*(1.+max(0., muon.MiniPFRelIso()-0.1));
        // For Run3, restrict the ptCorr < 50 GeV
        ptCorr = ((Run == 3) ? (ptCorr > 50. ? 49 : ptCorr) : ptCorr);
        const string cset_name = "fakerate_muon_" + string(syst_key.Data());
        auto cset = cset_muon_TopHNT_fakerate->at(cset_name);
        return safeEvaluate(cset, "GetFakeRate", {absEta, ptCorr});
    } else {
        // This is only a template
        auto cset = cset_muon->at(string(id_key));
        return safeEvaluate(cset, "GetFakeRate", {fabs(muon.Eta()), muon.OriginalPt()});
    }
}

// Electron
// For Run2, scale uncertainty is not stored in the NanoAODv9.
// Should patch from https://github.com/cms-egamma/ScaleFactorsJSON
// https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3
float MyCorrection::GetElectronScaleUnc(const float scEta, const unsigned char seedGain, const unsigned int runNumber, const float r9, const float pt, const variation syst) const {
    if (IsDATA) return 1.0;
    
    switch(Run) {
    case 2: {
        if (syst == variation::nom) return 1.;
        auto cset = cset_electron_variation->at("UL-EGM_ScaleUnc");
        vector<correction::Variable::Type> args = {
            GetEra().Data(),
            getSystString_EGMScale(syst),
            scEta,
            static_cast<int>(seedGain)
        };
        return safeEvaluate(cset, "GetElectronScaleSF", args);
    }
    case 3: {
        if (syst == variation::nom) return 1.;
        const string key = (GetEra().Contains("2022")) ? "Scale" : EGM_keys.at(GetEra().Data())+"_ScaleJSON";
        auto cset = cset_electron_variation->at(key);
        vector<correction::Variable::Type> args = {
            "total_uncertainty",
            static_cast<int>(seedGain),
            static_cast<float>(runNumber),
            scEta,
            r9,
            pt
        };
        const float unc = safeEvaluate(cset, "GetElectronScaleSF", args);
        if (syst == variation::up) return 1.+unc;
        else if (syst == variation::down) return 1.-unc;
        else
            throw runtime_error("[MyCorrection::GetElectronScaleUnc] Invalid syst value");
    }
    default:
        throw runtime_error("[MyCorrection::GetElectronScaleUnc] Invalid run number");
    }
    
    // This should never be reached, but added to avoid compiler warning
    return 1.0;
}
    
float MyCorrection::GetElectronSmearUnc(const Electron &electron, const variation syst, const unsigned int seed) const {
    if (IsDATA) return 1.0; // No smearing for data, only applied to MC
    if (Run == 2) throw runtime_error("[MyCorrection::GetElectronSmearUnc] Run2 is not supported by NanoAODv9");

    const string key = (GetEra().Contains("2022")) ? "Smearing" : EGM_keys.at(GetEra().Data())+"_SmearingJSON";
    auto cset = cset_electron_variation->at(key);
    vector<correction::Variable::Type> args = {
        "rho",
        electron.scEta(),
        electron.r9()
    };
    const float rho = safeEvaluate(cset, "GetElectronScaleSF", args);
    
    TRandom3 rng(seed);
    
    // Handle different variation cases
    if (syst == variation::nom) {
        // For nominal case, apply normal smearing
        return rng.Gaus(1.0, rho);
    } else if (syst == variation::up) {
        // For up variation, increase the width of the Gaussian
        return rng.Gaus(1.0, rho + safeEvaluate(cset, "GetElectronScaleSF", {"err_rho", electron.scEta(), electron.r9()}));
    } else if (syst == variation::down) {
        // For down variation, decrease the width of the Gaussian
        return rng.Gaus(1.0, rho - safeEvaluate(cset, "GetElectronScaleSF", {"err_rho", electron.scEta(), electron.r9()}));
    } else {
        throw runtime_error("[MyCorrection::GetElectronSmearUnc] Invalid syst value");
    }
}

float MyCorrection::GetElectronRECOSF(const float eta, const float pt, const float phi, const variation syst) const {
    switch(Run) {
        case 2:
            if (pt < 20.)
                return GetElectronIDSF("RecoBelow20", eta, pt, phi, syst);
            else
                return GetElectronIDSF("RecoAbove20", eta, pt, phi, syst);
            break;
        case 3:
            if (pt < 20.)
                return GetElectronIDSF("RecoBelow20", eta, pt, phi, syst);
            else if (pt < 75.)
                return GetElectronIDSF("Reco20to75", eta, pt, phi, syst);
            else
                return GetElectronIDSF("RecoAbove75", eta, pt, phi, syst);
            break;
        default:
            throw runtime_error("[MyCorrection::GetElectronRECOSF] Invalid run number");
    }
}

float MyCorrection::GetElectronRECOSF(const RVec<Electron> &electrons, const variation syst) const {
    float weight = 1.;
    for (const auto &electron : electrons) {
        weight *= GetElectronRECOSF(fabs(electron.Eta()), electron.Pt(), electron.Phi(), syst);
    }
    return weight;
}

float MyCorrection::GetElectronIDSF(const TString &Electron_ID_SF_Key, const float eta, const float pt, const float phi, const variation syst) const {
    if (Electron_ID_SF_Key == "TopHNT") {
        auto cset = cset_electron_TopHNT_idsf->at("sf");
        //const bool isLowStat = (GetEra() == "2022" || GetEra() == "2023BPix");
        const float maxPt = (Run == 3) ? 99. : 199.;
        if (syst == variation::nom) {
            return safeEvaluate(cset, "GetElectronIDSF", {eta, min(pt, maxPt), "nom"});
        } else if (syst == variation::up) {
            return safeEvaluate(cset, "GetElectronIDSF", {eta, min(pt, maxPt), "up"});
        } else if (syst == variation::down) {
            return safeEvaluate(cset, "GetElectronIDSF", {eta, min(pt, maxPt), "down"});
        } else {
            throw runtime_error("[MyCorrection::GetElectronIDSF] Invalid syst value");
        }
    } else {
        // POG IDs
        string key;
        if (Run == 2) key = "UL-Electron-ID-SF";
        else if (Run == 3) key = "Electron-ID-SF";
        else throw runtime_error("[MyCorrection::GetElectronIDSF] Invalid run number");
        
        auto cset = cset_electron->at(key);
        //NOTE: from 2023, It seems some SF depends on phi.
        if(!isInputInCorrection("phi", cset)) {
            return safeEvaluate(cset, "GetElectronIDSF", {EGM_keys.at(GetEra().Data()), getSystString_EGM(syst), string(Electron_ID_SF_Key), eta, pt});
        } else {
            return safeEvaluate(cset, "GetElectronIDSF", {EGM_keys.at(GetEra().Data()), getSystString_EGM(syst), string(Electron_ID_SF_Key), eta, pt, phi});
        }
    }
}

float MyCorrection::GetElectronIDSF(const TString &Electron_ID_SF_Key, const RVec<Electron> &electrons, const variation syst) const {
    float weight = 1.;
    for (const auto &electron : electrons) {
        if (Electron_ID_SF_Key == "TopHNT") {
            weight *= GetElectronIDSF(Electron_ID_SF_Key, electron.scEta(), electron.Pt(), 0., syst);
        } else {
            weight *= GetElectronIDSF(Electron_ID_SF_Key, fabs(electron.Eta()), electron.Pt(), electron.Phi(), syst);
        }
    }
    return weight;
}

float MyCorrection::GetFakeRate(const Electron &electron, const TString &id_key, const TString &syst_key) const {
    if (id_key == "TopHNT") {
        //const float scEta = (Run == 2) ? fabs(electron.scEta()): electron.scEta();
        const float absScEta = fabs(electron.scEta());
        const float ptCorr = electron.Pt()*(1.+max(0., electron.MiniPFRelIso()-0.1));
        const string cset_name = "fakerate_electron_" + string(syst_key.Data());
        auto cset = cset_electron_TopHNT_fakerate->at(cset_name);
        return safeEvaluate(cset, "GetFakeRate", {absScEta, ptCorr});
    } else {
        // This is only a template
        auto cset = cset_electron->at(string(id_key));
        return safeEvaluate(cset, "GetFakeRate", {fabs(electron.Eta()), electron.Pt()});
    }
}

// Trigger
float MyCorrection::GetMuonTriggerEff(const TString &Muon_Trigger_Eff_Key, const float abseta, const float pt, const bool isData, const variation syst) const {
    auto cset = cset_muon_trig_eff->at(string(Muon_Trigger_Eff_Key));
    if (isData)
        return safeEvaluate(cset, "GetTriggerEff", {"data", getSystString_MUO(syst), fabs(abseta), pt});
    else
        return safeEvaluate(cset, "GetTriggerEff", {"mc", getSystString_MUO(syst), fabs(abseta), pt});
}

float MyCorrection::GetMuonTriggerSF(const TString &Muon_Trigger_SF_Key, const RVec<Muon> &muons, const variation syst) const {
    float eff_data = 1.;
    float eff_mc = 1.;
    float weight = 1.;
    for (const auto &muon : muons) {
        eff_data *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()), muon.Pt(), true, syst);
        eff_mc *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()), muon.Pt(), false, syst);
    }
    weight = (1. - eff_data) / (1. - eff_mc);
    return weight;
}

float MyCorrection::GetElectronTriggerEff(const TString &Electron_Trigger_SF_Key, const float eta, const float pt, const float phi, const bool isDATA, const variation syst) const {
    string key = isDATA ? "Electron-HLT-DataEff" : "Electron-HLT-McEff";
    auto cset = cset_electron_hlt->at(key);
    // hardcoded replacemet
    string ValType = getSystString_EGM(syst);
    if (ValType == "sf")
        ValType = "nom";
    else if (ValType == "sfup")
        ValType = "up";
    else if (ValType == "sfdown")
        ValType = "down";
    else
        throw runtime_error("[MyCorrection::GetElectronTriggerEff] Invalid syst value");
    try {
        if(!isInputInCorrection("phi", cset)){
            return safeEvaluate(cset, "GetTriggerEff", {EGM_keys.at(GetEra().Data()), ValType, string(Electron_Trigger_SF_Key), eta, pt});
        } else {
            return safeEvaluate(cset, "GetTriggerEff", {EGM_keys.at(GetEra().Data()), ValType, string(Electron_Trigger_SF_Key), eta, pt, phi});
        }
    } catch (exception &e) {
        cerr << "[MyCorrection::GetElectronTriggerEff] " << e.what() << endl;
        throw e;
    }
}

float MyCorrection::GetElectronTriggerSF(const TString &Electron_Trigger_SF_Key, const float eta, const float pt, const float phi, const variation syst) const {
    auto cset = cset_electron_hlt->at("Electron-HLT-SF");
    try {
        if(!isInputInCorrection("phi", cset)){
            return safeEvaluate(cset, "GetTriggerEff", {EGM_keys.at(GetEra().Data()), getSystString_EGM(syst), string(Electron_Trigger_SF_Key), eta, pt});
        } else {
            return safeEvaluate(cset, "GetTriggerEff", {EGM_keys.at(GetEra().Data()), getSystString_EGM(syst), string(Electron_Trigger_SF_Key), eta, pt, phi});
        }
    } catch (exception &e) {
        cerr << "[MyCorrection::GetElectronTriggerSF] " << e.what() << endl;
        throw e;
    }
}

// double lepton triggers
// This function is used for leptons passing TopHNT ID
float MyCorrection::GetTriggerEff(const Muon &muon, const TString &trigger_leg_key, const bool isData, const variation syst) const {
    if (trigger_leg_key == "DblMu_Mu17Leg") {
        const string jsonkey = isData ? "data" : "sim";
        auto cset = cset_muon_TopHNT_dblmu_leg1_eff->at(jsonkey);
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(muon.Eta()), muon.OriginalPt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else if (trigger_leg_key == "DblMu_Mu8Leg") {
        const string jsonkey = isData ? "data" : "sim";
        auto cset = cset_muon_TopHNT_dblmu_leg2_eff->at(jsonkey);
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(muon.Eta()), muon.OriginalPt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else if (trigger_leg_key == "EMu_Mu23Leg") {
        const string jsonkey = isData ? "Mu23El12_Data" : "Mu23El12_MC";
        auto cset = cset_muon_TopHNT_emu_leg1_eff->at(jsonkey); 
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(muon.Eta()), muon.OriginalPt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else if (trigger_leg_key == "EMu_Mu8Leg") {
        const string jsonkey = isData ? "Mu8El23_Data" : "Mu8El23_MC";
        auto cset = cset_muon_TopHNT_emu_leg2_eff->at(jsonkey); 
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(muon.Eta()), muon.OriginalPt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else {
        throw runtime_error("[MyCorrection::GetTriggerEff] Invalid trigger leg key");
    }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetTriggerEff(const Electron &electron, const TString &trigger_leg_key, const bool isData, const variation syst) const {
    if (trigger_leg_key == "EMu_El23Leg") {
        const string jsonkey = isData ? "Mu8El23_Data" : "Mu8El23_MC";
        auto cset = cset_electron_TopHNT_emu_leg1_eff->at(jsonkey);
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(electron.scEta()), electron.Pt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else if (trigger_leg_key == "EMu_El12Leg") {
        const string jsonkey = isData ? "Mu23El12_Data" : "Mu23El12_MC";
        auto cset = cset_electron_TopHNT_emu_leg2_eff->at(jsonkey);
        float eff = safeEvaluate(cset, "GetTriggerEff", {fabs(electron.scEta()), electron.Pt(), getSystString_CUSTOM(syst)});
        return eff < 1. ? eff : 1.;
    } else {
        throw runtime_error("[MyCorrection::GetTriggerEff] Invalid trigger leg key");
    }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetDblMuTriggerEff(const RVec<Muon> &muons, const bool isData, const variation syst) const {
    const float filter_eff = GetPairwiseFilterEff("DblMu", isData);
    if (muons.size() == 2) {
        float leg1_eff = GetTriggerEff(muons.at(0), "DblMu_Mu17Leg", isData, syst);
        float leg2_eff = GetTriggerEff(muons.at(1), "DblMu_Mu8Leg", isData, syst);
        return leg1_eff*leg2_eff*filter_eff;
    } else if (muons.size() == 3) {
        const auto &mu1 = muons.at(0);
        const auto &mu2 = muons.at(1);
        const auto &mu3 = muons.at(2);

        float case1 = GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
              case1 *= GetTriggerEff(mu2, "DblMu_Mu8Leg", isData, syst);
              case1 *= filter_eff;
        float case2 = 1.-GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
              case2 *= GetTriggerEff(mu2, "DblMu_Mu17Leg", isData, syst);
              case2 *= GetTriggerEff(mu3, "DblMu_Mu8Leg", isData, syst);
              case2 *= filter_eff;
        float case3 = GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
              case3 *= (1. - GetTriggerEff(mu2, "DblMu_Mu8Leg", isData, syst))*filter_eff;
              case3 *= GetTriggerEff(mu3, "DblMu_Mu8Leg", isData, syst)*filter_eff;
        return case1 + case2 + case3;
    } else {
        throw runtime_error("[MyCorrection::GetDblMuTriggerEff] Invalid muon and electron size configuration");
    }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetDblMuTriggerSF(const RVec<Muon> &muons, const variation syst) const {
    float eff_data = -999.;
    float eff_mc = -999.;    
    if (syst == variation::nom) {
        eff_data = GetDblMuTriggerEff(muons, true, syst);
        eff_mc = GetDblMuTriggerEff(muons, false, syst);
    } else if (syst == variation::up) {
        eff_data = GetDblMuTriggerEff(muons, true, variation::up);
        eff_mc = GetDblMuTriggerEff(muons, false, variation::down);
    } else if (syst == variation::down) {
        eff_data = GetDblMuTriggerEff(muons, true, variation::down);
        eff_mc = GetDblMuTriggerEff(muons, false, variation::up);
    } else {
        throw runtime_error("[MyCorrection::GetDblMuTriggerSF] Invalid syst value");
    }
    return eff_mc > 0. ? eff_data / eff_mc : 1.;
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetEMuTriggerEff(const RVec<Electron> &electrons, const RVec<Muon> &muons, const bool isData, const variation syst) const {
    const float filter_eff = GetPairwiseFilterEff("EMu", isData);
    if (electrons.size() == 1 && muons.size() == 1) {
        const auto &el = electrons.at(0);
        const auto &mu = muons.at(0);
        float eff_el = mu.Pt() > 25. ? GetTriggerEff(el, "EMu_El12Leg", isData, syst) : GetTriggerEff(el, "EMu_El23Leg", isData, syst);
        float eff_mu = el.Pt() > 25. ? GetTriggerEff(mu, "EMu_Mu8Leg", isData, syst) : GetTriggerEff(mu, "EMu_Mu23Leg", isData, syst);
        return eff_el * eff_mu * filter_eff;
    } else if (electrons.size() == 1 && muons.size() == 2) {
        const auto &el = electrons.at(0);
        const auto &mu1 = muons.at(0);
        const auto &mu2 = muons.at(1);

        float case1 = GetTriggerEff(mu1, "EMu_Mu8Leg", isData, syst)
                    + (1.-GetTriggerEff(mu1, "EMu_Mu8Leg", isData, syst))*GetTriggerEff(mu2, "EMu_Mu8Leg", isData, syst)*filter_eff;
        float case2 = mu2.Pt() > 25. ? GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst)
                                        + (1.-GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst))*GetTriggerEff(mu2, "EMu_Mu23Leg", isData, syst)*filter_eff
                                    : GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst);
        float eff_el = (mu1.Pt() > 25. || mu2.Pt() > 25.) ? GetTriggerEff(el, "EMu_El12Leg", isData, syst) : GetTriggerEff(el, "EMu_El23Leg", isData, syst);
        float eff_mu = el.Pt() > 25. ? case1 : case2;
        return eff_el * eff_mu * filter_eff;
    } else {
        throw runtime_error("[MyCorrection::GetEMuTriggerEff] Invalid electron and muon size configuration");
    }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetEMuTriggerSF(const RVec<Electron> &electrons, const RVec<Muon> &muons, const variation syst) const {
    float eff_data = -999.;
    float eff_mc = -999.;    
    if (syst == variation::nom) {
        eff_data = GetEMuTriggerEff(electrons, muons, true, syst);
        eff_mc = GetEMuTriggerEff(electrons, muons, false, syst);
    } else if (syst == variation::up) {
        eff_data = GetEMuTriggerEff(electrons, muons, true, variation::up);
        eff_mc = GetEMuTriggerEff(electrons, muons, false, variation::down);
    } else if (syst == variation::down) {
        eff_data = GetEMuTriggerEff(electrons, muons, true, variation::down);
        eff_mc = GetEMuTriggerEff(electrons, muons, false, variation::up);
    } else {
        throw runtime_error("[MyCorrection::GetEMuTriggerSF] Invalid syst value");
    }
    return eff_mc > 0. ? eff_data / eff_mc : 1.;
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetPairwiseFilterEff(const TString &filter_name, const bool isData) const {
    if (filter_name.Contains("DblMu")) {
        if (GetEra() == "2016preVFP") {
            return 1.;
        } else if (GetEra() == "2016postVFP") {
            return isData? 0.9798 : 0.9968;
        } else if (GetEra() == "2017") {
            return isData? 0.9961 : 0.9958;
        } else if (GetEra() == "2018") {
            return isData? 0.9988 : 0.9998;
        } else if (GetEra() == "2022") {
            return isData? 0.9996 : 0.9998;
        } else if (GetEra() == "2022EE") {
            return isData? 0.9994 : 0.9998;
        } else if (GetEra() == "2023") {
            return isData? 0.9993 : 0.9997;
        } else if (GetEra() == "2023BPix") {
            return isData? 0.9989 : 0.9996;
        } else {
            cerr << "[MyCorrection::GetPairwiseFilterEff] " << filter_name << " is not implemented for " << GetEra() << endl;
            return 1.;
        }
    } else if (filter_name.Contains("EMu")) {
        if (GetEra() == "2016preVFP") {
            return 1.;
        } else if (GetEra() == "2016postVFP") {
            return isData? 0.9638 : 0.9878;
        } else if (GetEra() == "2017") {
            return isData? 0.9989 : 0.9955;
        } else if (GetEra() == "2018") {
            return isData? 0.9946 : 0.9981;
        } else if (GetEra() == "2022") {
            return isData? 0.9964 : 0.9985;
        } else if (GetEra() == "2022EE") {
            return isData? 0.9949 : 0.9976;
        } else if (GetEra() == "2023") {
            return isData? 0.9944 : 0.9977;
        } else if (GetEra() == "2023BPix") {
            return isData? 0.9930 : 0.9972;
        } else {
            cerr << "[MyCorrection::GetPairwiseFilterEff] " << filter_name << " is not implemented for " << GetEra() << endl;
            return 1.;
        }
    } else {
        throw runtime_error("[MyCorrection::GetPairwiseFilterEff] Invalid filter name");
    }
}


// Pileup
float MyCorrection::GetPUWeight(const float nTrueInt, const variation syst, const TString &source) const {
    // nota bene: Input should be nTrueInt, not nPileUp
    correction::Correction::Ref cset = nullptr;
    cset = cset_puWeights->at(LUM_keys.at(GetEra().Data()));
    try {
        return safeEvaluate(cset, "GetPUWeight", {nTrueInt, getSystString_LUM(syst)});
    } catch (exception &e) {
        cerr << "[MyCorrection::GetPUWeight] " << e.what() << endl;
        return 1.;
    }
}

// Heavy flavor tagging
void MyCorrection::SetTaggingParam(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) {
    global_tagger = tagger;
    global_wp = wp;
    global_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger);
    global_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp);
}

float MyCorrection::GetBTaggingWP() const {
    try {
        correction::Correction::Ref cset = cset_btagging->at(global_taggerStr + "_wp_values");
        return safeEvaluate(cset, "GetBTaggingWP", {global_wpStr});
    } catch (const exception &e) {
        cerr << "[Correction::GetBTaggingWP] Warning: Failed to evaluate WP '" << global_wpStr << "' for tagger '" << global_taggerStr << endl;
        throw runtime_error(e.what());
        return 1.f;
    }
}

float MyCorrection::GetBTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const {
    // Convert enumerations to strings
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr     = JetTagging::GetTaggerCorrectionWPStr(wp).Data();

    try {
        correction::Correction::Ref cset = cset_btagging->at(this_taggerStr + "_wp_values");
        return safeEvaluate(cset, "GetBTaggingWP", {this_wpStr});
    } catch (const exception &e) {
        cerr << "[Correction::GetBTaggingWP] Warning: Failed to evaluate WP '"
                  << this_wpStr << "' for tagger '" << this_taggerStr << endl;
        throw runtime_error(e.what());
        return 1.f;
    }
}

float MyCorrection::GetBTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst) {
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    auto cset = cset_btagging_eff->at(this_taggerStr);
    return safeEvaluate(cset, "GetBTaggingSF", {"central", this_wpStr, flav, fabs(eta), pt});
}

// For fixed point WP
// sources are central, hf_corr, hf_uncorr, lf_corr, lf_uncorr
float MyCorrection::GetBTaggingReweightMethod1a(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method, const variation syst, const TString &source) {
    if (tagger != JetTagging::JetFlavTagger::DeepJet) {
        cerr << "[MyCorrection::GetBTaggingSF] Current implementation is for DeepJet tagger << endl";
        exit(EXIT_FAILURE);
    }
    if (! (method == JetTagging::JetTaggingSFMethod::comb || method == JetTagging::JetTaggingSFMethod::mujets)) {
        cerr << "[MyCorrection::GetBTaggingSF] Current implementation is for comb or mujets method only" << endl;
        exit(EXIT_FAILURE);
    }

    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    string syst_str = getSystString_BTV(syst);
    
    // Set source string for hf / lf
    string source_lf = "central";
    string source_hf = "central";
    if (source == "hf_corr") {
        source_hf = syst_str+"_correlated";
    } else if (source == "hf_uncorr") {
        source_hf = syst_str+"_uncorrelated";
    } else if (source == "lf_corr") {
        source_lf = syst_str+"_correlated";
    } else if (source == "lf_uncorr") {
        source_lf = syst_str+"_uncorrelated";
    } else {
        try {
            assert(source == "central");
        } catch (const std::exception &e) {
            cerr << "[MyCorrection::GetBTaggingSF] Invalid source: " << source << endl;
            exit(EXIT_FAILURE);
        }
    }

    float weight = 1.;
    auto cset = cset_btagging->at(this_taggerStr+"_"+JetTagging::GetJetTaggingSFMethodStr(method).Data());
    string light_str = (Run == 2) ? this_taggerStr+"_incl" : this_taggerStr+"_light";
    auto cset_light = cset_btagging->at(light_str);
    float this_cut = GetBTaggingWP(tagger, wp);
    for (const auto &jet : jets) {
        const bool is_heavy = abs(jet.hadronFlavour()) == 5 || abs(jet.hadronFlavour()) == 4;
        const int this_flav = is_heavy ? abs(jet.hadronFlavour()) : 0;
        auto this_source = is_heavy ? source_hf : source_lf;
        auto this_cset = is_heavy ? cset : cset_light;
        const float eff = GetBTaggingEff(jet.Eta(), jet.Pt(), this_flav, tagger, wp, syst);

        const float sf = safeEvaluate(this_cset, "GetBTaggingSF", {this_source, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
        if (jet.GetBTaggerResult(tagger) > this_cut) {
            weight *= sf;
        } else {
            // To avoid dividing by 0
            const float denom = (eff == 1) ? 0.01 : (1.-eff);
            weight *= (1. - eff * sf) / denom;
        }
    }
    return weight;
}



float MyCorrection::GetBTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method, const variation syst, const TString &source) {
    if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet) {
        cerr << "[MyCorrection::GetBTaggingSF] DeepJet is the only supported tagger for Run2 UL" << endl;
        exit(EXIT_FAILURE);
    }
    
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    string syst_str = getSystString_BTV(syst);
    unordered_set<string> c_flav_source = {"cferr1", "cferr2"};
    if (source != "total")
        syst_str = getSystString_BTV(syst) + "_" + source;
    float weight = 1.;
    // method is comb, mujets, or shape
    if (method == JetTagging::JetTaggingSFMethod::shape) {
        auto cset = cset_btagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        for (const auto &jet : jets) {
            int this_flav = 0;
            if (abs(jet.hadronFlavour()) == 5)
                this_flav = 5;
            else if (abs(jet.hadronFlavour()) == 4)
                this_flav = 4;
            float this_score = jet.GetBTaggerResult(tagger);
            if (this_score < 0.f)
                continue; // defaulted jet
            if ((c_flav_source.find(source.Data()) != c_flav_source.end())) {
                weight *= safeEvaluate(cset, "GetBTaggingSF", {this_flav == 4 ? syst_str : getSystString_BTV(MyCorrection::variation::nom),
                                          this_flav,
                                          fabs(jet.Eta()),
                                          jet.Pt(),
                                          this_score});
            } else {
                weight *= safeEvaluate(cset, "GetBTaggingSF", {this_flav == 4 ? getSystString_BTV(MyCorrection::variation::nom) : syst_str,
                                          this_flav,
                                          fabs(jet.Eta()),
                                          jet.Pt(),
                                          this_score});
            }
        }
        return weight;
    } else if ((method == JetTagging::JetTaggingSFMethod::comb) || (method == JetTagging::JetTaggingSFMethod::mujets)) {
        auto cset = cset_btagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        string light_str;
        if (Run == 2)
            light_str = this_taggerStr + "_incl";
        else if (Run == 3)
            light_str = this_taggerStr + "_light";
        auto cset_light = cset_btagging->at(light_str);
        float this_cut = GetBTaggingWP(tagger, wp);
        for (const auto &jet : jets) {
            const bool is_heavy = abs(jet.hadronFlavour()) == 5 || abs(jet.hadronFlavour()) == 4;
            const int this_flav = is_heavy ? abs(jet.hadronFlavour()) : 0;

            const float eff = GetBTaggingEff(jet.Eta(), jet.Pt(), this_flav, tagger, wp, syst);
            auto this_cset = is_heavy ? cset : cset_light;

            const float sf = safeEvaluate(this_cset, "GetBTaggingSF", {syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
            if (jet.GetBTaggerResult(tagger) > this_cut) {
                weight *= sf;
            } else {
                weight *= (1. - eff * sf) / (1. - eff);
            }
        }
        return weight;
    } else {
        cout << "[MyCorrection::GetBTaggingSF] method " << JetTagging::GetJetTaggingSFMethodStr(method) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float MyCorrection::GetBTaggingR(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, string &processName, const variation syst, const TString &source) const {
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string syst_str = getSystString_BTV(syst);

    if (source != "total")
        syst_str = getSystString_BTV(syst) + "_" + source;
    float weight = 1.;

    this_taggerStr += (string("_") + processName);

    // First, check cset_btagging_R is loaded or nullptr
    if(!cset_btagging_R) throw std::runtime_error("[MyCorrection::GetBTaggingR] cset_btagging_R is not loaded or nullptr");
    auto cset = cset_btagging_R->at(this_taggerStr);

    for (const auto &jet : jets) {
        int this_flav = 0;
        if (abs(jet.hadronFlavour()) == 5)
            this_flav = 5;
        else if (abs(jet.hadronFlavour()) == 4)
            this_flav = 4;

        weight *= safeEvaluate(cset, "GetBTaggingEff", {syst_str,
                                  this_flav,
                                  jet.Pt(),
                                  fabs(jet.Eta())});
    }
    return weight;
}

pair<float, float> MyCorrection::GetCTaggingWP() const {
    try
    {
        correction::Correction::Ref cset = cset_ctagging->at(global_taggerStr + "_wp_values");
        float valCvB = safeEvaluate(cset, "GetCTaggingWP", {global_wpStr, "CvB"});
        float valCvL = safeEvaluate(cset, "GetCTaggingWP", {global_wpStr, "CvL"});
        return make_pair(valCvB, valCvL);
    }
    catch (const exception& e)
    {
        // If the requested WP is not found or any other error occurs, 
        // log a warning and return (1.f, 1.f) as a fallback.
        cerr << "[Correction::GetCTaggingWP] Warning: Failed to evaluate WP '"
             << global_wpStr << "' for tagger '" << global_taggerStr
             << endl;
        throw runtime_error(e.what());
        return make_pair(1.f, 1.f);
    }
}


pair<float, float> MyCorrection::GetCTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const
{
    // Convert enumerations to strings using your existing utility functions
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr     = JetTagging::GetTaggerCorrectionWPStr(wp).Data();

    // Retrieve the relevant correction set
    correction::Correction::Ref cset = cset_ctagging->at(this_taggerStr + "_wp_values");

    try
    {
        // Evaluate the corrections. If the WP does not exist, an exception might be thrown
        float valCvB = safeEvaluate(cset, "GetCTaggingWP", {this_wpStr, "CvB"});
        float valCvL = safeEvaluate(cset, "GetCTaggingWP", {this_wpStr, "CvL"});

        // If everything is fine, return the pair
        return make_pair(valCvB, valCvL);
    }
    catch (const exception& e)
    {
        // In case the WP is not found (or any other error occurs),
        // print a warning (optional) and return default values
        cerr << "[Correction::GetCTaggingWP] Warning: WP '" 
                  << this_wpStr << "' not found for tagger '"
                  << this_taggerStr << endl;
        throw runtime_error(e.what());
        return make_pair(1.f, 1.f);
    }
}


float MyCorrection::GetCTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst)
{
    return 1.;
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    correction::Correction::Ref cset = cset_btagging_eff->at(this_taggerStr);
    return safeEvaluate(cset, "GetBTaggingEff", {
        getSystString_BTV(syst),
        this_wpStr,
    });
}

float MyCorrection::GetCTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod method, const variation syst, const TString &source)
{
    if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet)
    {
        cerr << "[MyCorrection::GetCTaggingSF] Run2 only supports DeepJet tagger" << endl;
        return 1.;
    }
    if (Run == 3)
    {
        // cerr << "[MyCorrection::GetCTaggingSF] Run3 C-Tagger SF not provided by POG yet" << endl;
        return 1.;
    }
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    string syst_str = getSystString_BTV(syst);
    if (source != "total")
        syst_str = getSystString_BTV(syst) + "_" + source;
    float weight = 1.;
    // method is comb, mujets, or shape
    if (method == JetTagging::JetTaggingSFMethod::shape)
    {
        auto cset = cset_ctagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        for (const auto &jet : jets)
        {
            int this_flav = 0;
            if (abs(jet.hadronFlavour()) == 5)
                this_flav = 5;
            else if (abs(jet.hadronFlavour()) == 4)
                this_flav = 4;
            pair<float, float> this_ctag = jet.GetCTaggerResult(tagger);
            if (this_ctag.first < 0.f || this_ctag.second < 0.f)
                continue; // defaulted jet
            weight *= safeEvaluate(cset, "GetCTaggingSF", {syst_str,
                                      this_flav,
                                      this_ctag.second,
                                      this_ctag.first});
        }
        return weight;
    }
    else if (method == JetTagging::JetTaggingSFMethod::wp)
    {
        auto cset = cset_btagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        auto cset_light = cset_btagging->at(this_taggerStr + "_incl");
        pair<float, float> this_cut = GetCTaggingWP(tagger, wp);
        pair<float, float> this_score;
        for (const auto &jet : jets)
        {
            this_score = jet.GetCTaggerResult(tagger);
            int this_flav = 0;
            if (abs(jet.hadronFlavour()) == 5)
                this_flav = 5;
            else if (abs(jet.hadronFlavour()) == 4)
                this_flav = 4;
            float eff = GetCTaggingEff(jet.Eta(), jet.Pt(), this_flav, tagger, wp, syst);
            float sf;
            if (this_flav == 0)
                sf = cset_light->evaluate({syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
            else
                sf = safeEvaluate(cset, "GetCTaggingSF", {syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
            if (this_score.first > this_cut.first && this_score.second > this_cut.second)
            {
                weight *= sf;
            }
            else
            {
                weight *= (1. - eff * sf) / (1. - eff);
            }
        }
        return weight;
    }
    else
    {
        cout << "[MyCorrection::GetCTaggingSF] method " << JetTagging::GetJetTaggingSFMethodStr(method) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float MyCorrection::GetCTaggingR(const float nTrueInt, const float HT, const JetTagging::JetFlavTagger tagger, const TString &processName, const TString &ttBarCategory, const TString &syst_str) const
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    if (processName != "")
        this_taggerStr += "_" + processName;
    else
        this_taggerStr += (string("_") + Sample.Data());
    auto cset = cset_ctagging_R->at(this_taggerStr);
    return safeEvaluate(cset, "GetTopPtReweight", {syst_str.Data(), ttBarCategory.Data(), nTrueInt, HT});
}

// Pileup Jet ID
float MyCorrection::GetPileupJetIDSF(const RVec<Jet> &jets, const unordered_map<int, int> &matched_idx, const TString &wp, const variation syst) {
    // Should pass jets after PUID, no mistag rate correction
    if (Run == 3) return 1.;

    float weight = 1.;
    auto cset = cset_jmar->at(JME_PILEUP_keys.at(GetEra().Data()));
    string wp_str;
    if (wp == "tight") wp_str = "T";
    else if (wp == "medium") wp_str = "M";
    else if (wp == "loose") wp_str = "L";
    else throw std::invalid_argument("[MCCorrection::GetPileupJetIDSF] Invalid WP");

    for (int i = 0; i < jets.size(); i++) {
        if (jets.at(i).Pt() > 50.) continue;
        if (matched_idx.find(i) == matched_idx.end() || matched_idx.at(i) < 0) continue; // not matched
        float this_eff = safeEvaluate(cset, "GetJetVetoMapEff", {jets.at(i).Eta(), jets.at(i).Pt(), getSystString_JME(syst), wp_str});
        weight *= this_eff;
    }
    return weight;
}

// JERC
float MyCorrection::GetJER(const float eta, const float pt, const float rho) const {
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JER_GT.at(GetEra().Data());
    cset_string.replace(cset_string.find("######"), 6, "PtResolution");
    cset = cset_jerc->at(cset_string);
    return safeEvaluate(cset, "GetJER", {eta, pt, rho});
}

float MyCorrection::GetJERSF(const float eta, const float pt, const variation syst, const TString &source) const {
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JER_GT.at(GetEra().Data());
    cset_string.replace(cset_string.find("######"), 6, "ScaleFactor");
    cset = cset_jerc->at(cset_string);
    if (Run == 3)
    {
        return safeEvaluate(cset, "GetJERSF", {eta, pt, getSystString_JME(syst)});
    }
    else if (Run == 2)
    {
        return safeEvaluate(cset, "GetJERSF", {eta, getSystString_JME(syst)});
    }
    return 1.;
}

//JESC
float MyCorrection::GetJESSF(const float area, const float eta, const float pt, const float phi, const float rho, const unsigned int runNumber) const {
    correction::CompoundCorrection::Ref cset = nullptr;
    string cset_string = JME_JES_GT.at(GetEra().Data());
    cset_string.replace(cset_string.find("######"), 6, "L1L2L3Res");
    cset = cset_jerc->compound().at(cset_string);
    vector<correction::Variable::Type> args;
    float JESSF = 1.;
    if (GetEra() == "2023BPix" || GetEra() == "2024") {
        args = {area, eta, pt, rho, phi};
        if (IsDATA) args = {area, eta, pt, rho, phi, static_cast<float>(runNumber)};
    } else if (GetEra() == "2023") {
        args = {area, eta, pt, rho};
        if (IsDATA) args = {area, eta, pt, rho, static_cast<float>(runNumber)};
    } else {
        args = {area, eta, pt, rho};
    }
    return safeEvaluate(cset, "GetJERSF", args);
}

float MyCorrection::GetJESUncertainty(const float eta, const float pt, const variation syst, const TString &source) const {
    int int_syst = 0;
    if (syst == variation::up)
        int_syst = 1;
    else if (syst == variation::down)
        int_syst = -1;
    else
        int_syst = 0;

    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JES_GT.at(GetEra().Data());
    cset_string.replace(cset_string.find("######"), 6, source);
    cset = cset_jerc->at(cset_string);
    float this_factor = 1.;
    this_factor += (int_syst * safeEvaluate(cset, "GetJESUncertainty", {eta, pt}));
    return this_factor;
}

bool MyCorrection::IsJetVetoZone(const float eta, const float phi, TString mapCategory) const
{
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_vetomap_keys.at(GetEra().Data());
    cset = cset_jetvetomap->at(cset_string);
    if (safeEvaluate(cset, "IsJetVetoZone", {mapCategory.Data(), eta, phi}) > 0) return true;
    return false;
}

void MyCorrection::METXYCorrection(Particle &Met, const int RunNumber, const int npvs, const XYCorrection_MetType MetType) {
    if (Run == 3) {
        // Not supported yet
        return;
    }

    correction::Correction::Ref cset_pt = nullptr;
    correction::Correction::Ref cset_phi = nullptr;
    switch (MetType) {
    case XYCorrection_MetType::Type1PFMET:
        if (IsDATA) {
            cset_pt = cset_met->at("pt_metphicorr_pfmet_data");
            cset_phi = cset_met->at("phi_metphicorr_pfmet_data");
        } else {
            cset_pt = cset_met->at("pt_metphicorr_pfmet_mc");
            cset_phi = cset_met->at("phi_metphicorr_pfmet_mc");
        }
        break;
    case XYCorrection_MetType::Type1PuppiMET:
        cout << "[MyCorrection::METXYCorrection] Type1PuppiMET is not recommended for PUPPI MET" << endl;
        if (IsDATA) {
            cset_pt = cset_met->at("pt_metphicorr_puppimet_data");
            cset_phi = cset_met->at("phi_metphicorr_puppimet_data");
        } else {
            cset_pt = cset_met->at("pt_metphicorr_puppimet_mc");
            cset_phi = cset_met->at("phi_metphicorr_puppimet_mc");
        }
        break;
    }
    float this_pt = cset_pt->evaluate({Met.Pt(), Met.Phi(), static_cast<float>(npvs), static_cast<float>(RunNumber)});
    float this_phi = cset_phi->evaluate({Met.Pt(), Met.Phi(), static_cast<float>(npvs), static_cast<float>(RunNumber)});
    float this_eta = Met.Eta();
    float this_m = Met.M();
    Met.SetPtEtaPhiM(this_pt, this_eta, this_phi, this_m);
}


float MyCorrection::GetTopPtReweight(const RVec<Gen> &gens) const {
    if (!Sample.Contains("TT") || !Sample.Contains("powheg"))
        return 1.;
    const Gen *top = nullptr;
    const Gen *antitop = nullptr;
    // loop while find top and antitop
    for (const auto &gen : gens)
    {
        if (gen.PID() == 6 and gen.isHardProcess())
            top = &gen;
        if (gen.PID() == -6 and gen.isHardProcess())
            antitop = &gen;
        if (top != nullptr && antitop != nullptr)
            break;
    }
    if (top == nullptr || antitop == nullptr)
        return 1.;

    if (top->Pt() > 2000 || antitop->Pt() > 2000)
        return 1.;

    float a = 0.103;
    float b = -0.0118;
    float c = 0.000134;
    float d = 0.973;
    float top_sf = a * TMath::Exp(b * top->Pt()) - c * top->Pt() + d;
    float antitop_sf = a * TMath::Exp(b * antitop->Pt()) - c * antitop->Pt() + d;
    return TMath::Sqrt(top_sf * antitop_sf);
}

float MyCorrection::GetWZNjetsSF(const float &njets, const string &syst) const {
    if (! (Sample.Contains("WZTo3LNu"))) return 1.;

    string cset_string;
    if (Run == 2) {
        cset_string = "WZNjetsSF_WZCombined_Run2_" + syst;
    } else if (Run == 3) {
        cset_string = "WZNjetsSF_WZCombined_Run3_" + syst;
    } else {
        throw invalid_argument("[MyCorrection::GetWZNjetsSF] Invalid run number");
    }
    auto cset = cset_WZ_Nj_SF->at(cset_string);
    return safeEvaluate(cset, "GetWZNjetsSF", {njets});
}
