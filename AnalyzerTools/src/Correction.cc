#include "Correction.h"

Correction::Correction() {}

Correction::Correction(const TString &era, const TString &sample, const bool IsData, const string &btagging_eff_file, const string &ctagging_eff_file, const string &btagging_R_file, const string &ctagging_R_file)
{
    cout << "[Correction::Correction] MCorrection created for " << era << endl;
    SetEra(era);
    SetSample(sample);
    setIsData(IsData);
    string json_muon = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/MUO";
    string json_muon_trig_eff = string(getenv("SKNANO_DATA"));
    string json_puWeights = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/LUM";
    string json_btagging = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/BTV";
    string json_ctagging = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/BTV";
    string json_btagging_eff = string(getenv("SKNANO_DATA"));
    string json_ctagging_eff = string(getenv("SKNANO_DATA"));
    string json_btagging_R = string(getenv("SKNANO_DATA"));
    string json_ctagging_R = string(getenv("SKNANO_DATA"));
    string json_electron = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/EGM";
    string json_electron_hlt = json_electron;
    string json_photon = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/EGM";
    string json_jerc = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/JME";
    string json_jerc_fatjet = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/JME";
    string json_jetvetomap = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/JME";
    string json_met = string(getenv("JSONPOG_INTEGRATION_PATH")) + "/POG/JME";

    if (era == "2016preVFP")
    {
        json_muon += "/2016preVFP_UL/muon_Z.json.gz";
        json_muon_trig_eff += "/2016preVFP/MUO/muon_trig.json";
        json_puWeights += "/2016preVFP_UL/puWeights.json.gz";
        json_btagging += "/2016preVFP_UL/btagging.json.gz";
        json_ctagging += "/2016preVFP_UL/ctagging.json.gz";
        json_btagging_eff += "/2016preVFP/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2016preVFP/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2016preVFP/BTV/" + btagging_R_file;
        json_ctagging_R += "/2016preVFP/BTV/" + ctagging_R_file;
        json_electron += "/2016preVFP_UL/electron.json.gz";
        json_photon += "/2016preVFP_UL/photon.json.gz";
        json_jerc += "/2016preVFP_UL/jet_jerc.json.gz";
        json_jerc_fatjet += "/2016preVFP_UL/fatJet_jerc.json.gz";
        json_jetvetomap += "/2016preVFP_UL/jetvetomaps.json.gz";
        json_met += "/2016preVFP_UL/met.json.gz";
    }
    else if (era == "2016postVFP")
    {
        json_muon += "/2016postVFP_UL/muon_Z.json.gz";
        json_muon_trig_eff += "/2016postVFP/MUO/muon_trig.json";
        json_puWeights += "/2016postVFP_UL/puWeights.json.gz";
        json_btagging += "/2016postVFP_UL/btagging.json.gz";
        json_ctagging += "/2016postVFP_UL/ctagging.json.gz";
        json_btagging_eff += "/2016postVFP/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2016postVFP/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2016postVFP/BTV/" + btagging_R_file;
        json_ctagging_R += "/2016postVFP/BTV/" + ctagging_R_file;
        json_electron += "/2016postVFP_UL/electron.json.gz";
        json_photon += "/2016postVFP_UL/photon.json.gz";
        json_jerc += "/2016postVFP_UL/jet_jerc.json.gz";
        json_jerc_fatjet += "/2016postVFP_UL/fatJet_jerc.json.gz";
        json_jetvetomap += "/2016postVFP_UL/jetvetomaps.json.gz";
        json_met += "/2016postVFP_UL/met.json.gz";
    }
    else if (era == "2017")
    {
        json_muon += "/2017_UL/muon_Z.json.gz";
        json_muon_trig_eff += "/2017/MUO/muon_trig.json";
        json_puWeights += "/2017_UL/puWeights.json.gz";
        json_btagging += "/2017_UL/btagging.json.gz";
        json_ctagging += "/2017_UL/ctagging.json.gz";
        json_btagging_eff += "/2017/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2017/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2017/BTV/" + btagging_R_file;
        json_ctagging_R += "/2017/BTV/" + ctagging_R_file;
        json_electron += "/2017_UL/electron.json.gz";
        json_photon += "/2017_UL/photon.json.gz";
        json_jerc += "/2017_UL/jet_jerc.json.gz";
        json_jerc_fatjet += "/2017_UL/fatJet_jerc.json.gz";
        json_jetvetomap += "/2017_UL/jetvetomaps.json.gz";
        json_met += "/2017_UL/met.json.gz";
    }
    else if (era == "2018")
    {
        json_muon += "/2018_UL/muon_Z.json.gz";
        json_muon_trig_eff += "/2018/MUO/muon_trig.json";
        json_puWeights += "/2018_UL/puWeights.json.gz";
        json_btagging += "/2018_UL/btagging.json.gz";
        json_ctagging += "/2018_UL/ctagging.json.gz";
        json_btagging_eff += "/2018/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2018/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2018/BTV/" + btagging_R_file;
        json_ctagging_R += "/2018/BTV/" + ctagging_R_file;
        json_electron += "/2018_UL/electron.json.gz";
        json_photon += "/2018_UL/photon.json.gz";
        json_jerc += "/2018_UL/jet_jerc.json.gz";
        json_jerc_fatjet += "/2018_UL/fatJet_jerc.json.gz";
        json_jetvetomap += "/2018_UL/jetvetomaps.json.gz";
        json_met += "/2018_UL/met.json.gz";
    }
    else if (era == "2022")
    {
        json_muon += "/2022_27Jun2023/muon_Z.json.gz";
        json_muon_trig_eff += "/2022/MUO/muon_trig.json";
        json_puWeights += "/2022_Summer22/puWeights.json.gz";
        json_btagging += "/2022_Summer22/btagging.json.gz";
        json_ctagging += "/2022_Summer22/ctagging.json.gz";
        json_btagging_eff += "/2022/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2022/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2022/BTV/" + btagging_R_file;
        json_ctagging_R += "/2022/BTV/" + ctagging_R_file;
        json_electron += "/2022_Summer22/electron.json.gz";
        json_electron_hlt += "/2022_Summer22/electronHlt.json.gz";
        json_photon += "/2022_Summer22/photon.json.gz";
        json_jerc += "/2022_Summer22/jet_jerc.json.gz";
        json_jerc_fatjet += "/2022_Summer22/fatJet_jerc.json.gz";
        json_jetvetomap += "/2022_Summer22/jetvetomaps.json.gz";
        json_met += "/2022_Summer22/met.json.gz";
    }
    else if (era == "2022EE")
    {
        json_muon += "/2022EE_27Jun2023/muon_Z.json.gz";
        json_muon_trig_eff += "/2022EE_27Jun2023/MUO/muon_trig.json";
        json_puWeights += "/2022_Summer22EE/puWeights.json.gz";
        json_btagging += "/2022_Summer22EE/btagging.json.gz";
        json_ctagging += "/2022_Summer22EE/ctagging.json.gz";
        json_btagging_eff += "/2022EE/BTV/" + btagging_eff_file;
        json_ctagging_eff += "/2022EE/BTV/" + ctagging_eff_file;
        json_btagging_R += "/2022EE/BTV/" + btagging_R_file;
        json_ctagging_R += "/2022EE/BTV/" + ctagging_R_file;
        json_electron += "/2022_Summer22EE/electron.json.gz";
        json_electron_hlt += "/2022_Summer22EE/electronHlt.json.gz";
        json_photon += "/2022_Summer22EE/photon.json.gz";
        json_jerc += "/2022_Summer22EE/jet_jerc.json.gz";
        json_jerc_fatjet += "/2022_Summer22EE/fatJet_jerc.json.gz";
        json_jetvetomap += "/2022_Summer22EE/jetvetomaps.json.gz";
        json_met += "/2022_Summer22EE/met.json.gz";
    }
    else
    {
        cerr << "[Correction::Correction] no correction file for era " << era << endl;
    }
    cout << "[Correction::Correction] using muon SF file: " << json_muon << endl;
    cset_muon = CorrectionSet::from_file(json_muon);

    cout << "[Correction::Correction] using puWeights file: " << json_puWeights << endl;
    cset_puWeights = CorrectionSet::from_file(json_puWeights);
    cout << "[Correction::Correction] using btagging file: " << json_btagging << endl;
    cset_btagging = CorrectionSet::from_file(json_btagging);
    cout << "[Correction::Correction] using ctagging file: " << json_ctagging << endl;
    cset_ctagging = CorrectionSet::from_file(json_ctagging);
    cout << "[Correction::Correction] using btagging eff file: " << json_btagging_eff << endl;
    cset_btagging_eff = CorrectionSet::from_file(json_btagging_eff);
    cout << "[Correction::Correction] using ctagging eff file: " << json_ctagging_eff << endl;
    cset_ctagging_eff = CorrectionSet::from_file(json_ctagging_eff);
    cout << "[Correction::Correction] using electron file: " << json_electron << endl;
    cset_electron = CorrectionSet::from_file(json_electron);
    cout << "[Correction::Correction] using photon file: " << json_photon << endl;
    cset_photon = CorrectionSet::from_file(json_photon);
    cout << "[Correction::Correction] using jerc file: " << json_jerc << endl;
    cset_jerc = CorrectionSet::from_file(json_jerc);
    cout << "[Correction::Correction] using jerc_fatjet file: " << json_jerc_fatjet << endl;
    cset_jerc_fatjet = CorrectionSet::from_file(json_jerc_fatjet);
    cout << "[Correction::Correction] using jetvetomap file: " << json_jetvetomap << endl;
    cset_jetvetomap = CorrectionSet::from_file(json_jetvetomap);
    // try open met file. if it does not exist, just skip it.
    try
    {
        cout << "[Correction::Correction] using muon trig eff file: " << json_muon_trig_eff << endl;
        cset_muon_trig_eff = CorrectionSet::from_file(json_muon_trig_eff);
    }
    catch (const std::exception &e)
    {
        cerr << "[Correction::Correction] Error with muon trig eff file: " << e.what() << endl;
    }

    try
    {
        cout << "[Correction::Correction] using electron hlt file: " << json_electron_hlt << endl;
        cset_electron_hlt = CorrectionSet::from_file(json_electron_hlt);
    }
    catch (const std::exception &e)
    {
        cerr << "[Correction::Correction] Error with electron hlt file: " << e.what() << endl;
    }

    try
    {
        cout << "[Correction::Correction] using met file: " << json_met << endl;
        cset_met = CorrectionSet::from_file(json_met);
    }
    catch (const std::exception &e)
    {
        cerr << "[Correction::Correction] Error with met file: " << e.what() << endl;
    }

    try
    {
        cout << "[Correction::Correction] using btagging R file: " << json_btagging_R << endl;
        cset_btagging_R = CorrectionSet::from_file(json_btagging_R);
    }
    catch (const std::exception &e)
    {
        cerr << "[Correction::Correction] Error with btagging R file: " << e.what() << endl;
    }

    try
    {
        cout << "[Correction::Correction] using ctagging R file: " << json_ctagging_R << endl;
        cset_ctagging_R = CorrectionSet::from_file(json_ctagging_R);
    }
    catch (const std::exception &e)
    {
        cerr << "[Correction::Correction] Error with ctagging R file: " << e.what() << endl;
    }

    LUM_keys["2022EE"] = "Collisions2022_359022_362760_eraEFG_GoldenJson";
    LUM_keys["2022"] = "Collisions2022_355100_357900_eraBCD_GoldenJson";
    LUM_keys["2018"] = "Collisions18_UltraLegacy_goldenJSON";
    LUM_keys["2017"] = "Collisions17_UltraLegacy_goldenJSON";
    LUM_keys["2016postVFP"] = "Collisions16_UltraLegacy_goldenJSON";
    LUM_keys["2016preVFP"] = "Collisions16_UltraLegacy_goldenJSON";

    EGM_keys["2022EE"] = "2022Re-recoE+PromptFG";
    EGM_keys["2022"] = "2022ReRecoBCD";
    EGM_keys["2018"] = "UL-Electron-ID-SF";
    EGM_keys["2017"] = "UL-Electron-ID-SF";
    EGM_keys["2016postVFP"] = "UL-Electron-ID-SF";
    EGM_keys["2016preVFP"] = "UL-Electron-ID-SF";

    JME_JER_GT["2022"] = "Summer22_22Sep2023_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2022EE"] = "Summer22EE_22Sep2023_JRV1_MC_######_AK4PFPuppi";
    JME_JER_GT["2018"] = "Summer19UL18_JRV2_MC_######_AK4PFchs";
    JME_JER_GT["2017"] = "Summer19UL17_JRV2_MC_######_AK4PFchs";
    JME_JER_GT["2016postVFP"] = "Summer20UL16_JRV3_MC_######_AK4PFchs";
    JME_JER_GT["2016preVFP"] = "Summer20UL16APV_JRV3_MC_######_AK4PFchs";

    JME_JES_GT["2022"] = "Summer22_22Sep2023_V2_MC_######_AK4PFPuppi";
    JME_JES_GT["2022EE"] = "Summer22EE_22Sep2023_V2_MC_######_AK4PFPuppi";
    JME_JES_GT["2018"] = "Summer19UL18_V5_MC_######_AK4PFchs";
    JME_JES_GT["2017"] = "Summer19UL17_V5_MC_######_AK4PFchs";
    JME_JES_GT["2016postVFP"] = "Summer19UL16APV_V7_MC_######_AK4PFchs";
    JME_JES_GT["2016preVFP"] = "Summer19UL16_V7_MC_V5_MC_######_AK4PFchs";

    JME_vetomap_keys["2022EE"] = "Summer22EE_23Sep2023_RunEFG_V1";
    JME_vetomap_keys["2022"] = "Summer22_23Sep2023_RunBCD_V1";
    JME_vetomap_keys["2018"] = "Summer19UL18_V1";
    JME_vetomap_keys["2017"] = "Summer19UL17_V1";
    JME_vetomap_keys["2016postVFP"] = "Summer19UL16_V1";
    JME_vetomap_keys["2016preVFP"] = "Summer19UL16_V1";
}

Correction::~Correction() {}

float Correction::GetMuonIDSF(const TString &Muon_ID_SF_Key, const float eta, const float pt, const variation syst, const TString &source) const
{
    auto cset = cset_muon->at(string(Muon_ID_SF_Key));
    return cset->evaluate({fabs(eta), pt, getSystString_MUO(syst)});
}

float Correction::GetMuonIDSF(const TString &Muon_ID_SF_Key, const RVec<Muon> &muons, const variation syst, const TString &source) const
{
    float weight = 1.;
    for (const auto &muon : muons)
    {
        weight *= GetMuonIDSF(Muon_ID_SF_Key, fabs(muon.Eta()), muon.Pt(), syst);
    }
    return weight;
}

float Correction::GetMuonTriggerEff(const TString &Muon_Trigger_Eff_Key, const float abseta, const float pt, const bool eff_for_data, const variation syst, const TString &source) const
{
    auto cset = cset_muon_trig_eff->at(string(Muon_Trigger_Eff_Key));
    if (eff_for_data)
        return cset->evaluate({"data", getSystString_MUO(syst), fabs(abseta), pt});
    else
        return cset->evaluate({"mc", getSystString_MUO(syst), fabs(abseta), pt});
}

float Correction::GetMuonTriggerWeight(const TString &Muon_Trigger_SF_Key, const RVec<Muon> &muons, const variation syst, const TString &source) const
{
    float eff_data = 1.;
    float eff_mc = 1.;
    float weight = 1.;
    for (const auto &muon : muons)
    {
        eff_data *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()), muon.Pt(), true, syst);
        eff_mc *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()), muon.Pt(), false, syst);
    }
    weight = (1. - eff_data) / (1. - eff_mc);
    return weight;
}

float Correction::GetElectronIDSF(const TString &Electron_ID_SF_Key, const float eta, const float pt, const variation syst, const TString &source) const
{
    auto cset = cset_electron->at("Electron-ID-SF");
    try
    {
        return cset->evaluate({EGM_keys.at(DataEra.Data()), getSystString_EGM(syst), string(Electron_ID_SF_Key), eta, pt});
    }
    catch (exception &e)
    {
        cerr << "[Correction::GetElectronIDSF] " << e.what() << endl;
        throw e;
    }
}

float Correction::GetElectronIDSF(const TString &Electron_ID_SF_Key, const RVec<Electron> &electrons, const variation syst, const TString &source) const
{
    float weight = 1.;
    for (const auto &electron : electrons)
    {
        weight *= GetElectronIDSF(Electron_ID_SF_Key, fabs(electron.Eta()), electron.Pt(), syst, source);
    }
    return weight;
}

float Correction::GetElectronTriggerSF(const TString &Electron_Trigger_SF_Key, const float eta, const float pt, const variation syst, const TString &source) const
{
    auto cset = cset_electron_hlt->at("Electron-HLT-SF");
    try
    {
        return cset->evaluate({EGM_keys.at(DataEra.Data()), getSystString_EGM(syst), string(Electron_Trigger_SF_Key), eta, pt});
    }
    catch (exception &e)
    {
        cerr << "[Correction::GetElectronTriggerSF] " << e.what() << endl;
        throw e;
    }
}

float Correction::GetElectronTriggerEff(const TString &Electron_Trigger_SF_Key, const float eta, const float pt, bool ofDATA, const variation syst, const TString &source) const
{
    std::string key = ofDATA ? "Electron-HLT-DataEff" : "Electron-HLT-McEff";
    auto cset = cset_electron_hlt->at(key);
    //hardcoded replacemet
    std::string ValType = getSystString_EGM(syst);
    if(ValType == "sf") ValType = "nom";
    else if(ValType == "sfup") ValType = "up";
    else if(ValType == "sfdown") ValType = "down";
    else throw std::runtime_error("[Correction::GetElectronTriggerEff] Invalid syst value");
    try
    {
        return cset->evaluate({EGM_keys.at(DataEra.Data()), ValType, string(Electron_Trigger_SF_Key), eta, pt});
    }
    catch (exception &e)
    {
        cerr << "[Correction::GetElectronTriggerEff] " << e.what() << endl;
        throw e;
    }
}

float Correction::GetElectronRECOSF(const float eta, const float pt, const variation syst, const TString &source) const
{
    if (pt < 20.)
        return GetElectronIDSF("RecoBelow20", eta, pt, syst);
    else if (pt >= 20. && pt < 75)
        return GetElectronIDSF("Reco20to75", eta, pt, syst);
    else if (pt >= 75.)
        return GetElectronIDSF("RecoAbove75", eta, pt, syst);
    else
        return 1.;
}

float Correction::GetElectronRECOSF(const RVec<Electron> &electrons, const variation syst, const TString &source) const
{
    float weight = 1.;
    for (const auto &electron : electrons)
    {
        weight *= GetElectronRECOSF(fabs(electron.Eta()), electron.Pt(), syst);
    }
    return weight;
}

float Correction::GetPUWeight(const float nTrueInt, const variation syst, const TString &source) const
{
    // nota bene: Input should be nTrueInt, not nPileUp
    correction::Correction::Ref cset = nullptr;
    cset = cset_puWeights->at(LUM_keys.at(DataEra.Data()));
    try
    {
        return cset->evaluate({nTrueInt, getSystString_LUM(syst)});
    }
    catch (exception &e)
    {
        cerr << "[Correction::GetPUWeight] " << e.what() << endl;
        return 1.;
    }
}
void Correction::SetTaggingParam(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp)
{
    global_tagger = tagger;
    global_wp = wp;
    global_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger);
    global_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp);
}

float Correction::GetBTaggingWP() const
{
    correction::Correction::Ref cset = cset_btagging->at(global_taggerStr + "_wp_values");
    return cset->evaluate({global_wpStr});
}

float Correction::GetBTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    correction::Correction::Ref cset = cset_btagging->at(this_taggerStr + "_wp_values");
    return cset->evaluate({this_wpStr});
}

float Correction::GetBTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst)
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    correction::Correction::Ref cset = cset_btagging_eff->at(this_taggerStr);
    return cset->evaluate({"eff", this_wpStr, flav, fabs(eta), pt});
}

float Correction::GetBTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod &method, const variation syst, const TString &source)
{
    if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet)
    {
        cerr << "[Correction::GetBTaggingSF] DeepJet is the only supported tagger for 2016preVFP, 2016postVFP, 2017, 2018, and 2018UL" << endl;
        return 1.;
    }
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    string syst_str = getSystString_BTV(syst);
    unordered_set<std::string> c_flav_source = {"cferr1",
                                           "cferr2"};
    if (source != "total")
        syst_str = getSystString_BTV(syst) + "_" + source;
    float weight = 1.;
    // method is comb, mujets, or shape
    if (method == JetTagging::JetTaggingSFMethod::shape)
    {
        auto cset = cset_btagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        for (const auto &jet : jets)
        {
            int this_flav = 0;
            if (abs(jet.hadronFlavour()) == 5)
                this_flav = 5;
            else if (abs(jet.hadronFlavour()) == 4)
                this_flav = 4;
            float this_score = jet.GetBTaggerResult(tagger);
            if (this_score < 0.f)
                continue; // defaulted jet
            if ((c_flav_source.find(source.Data()) != c_flav_source.end()))
            {
                weight *= cset->evaluate({this_flav == 4 ? syst_str : getSystString_BTV(Correction::variation::nom),
                                          this_flav,
                                          fabs(jet.Eta()),
                                          jet.Pt(),
                                          this_score});
            }
            else
            {
                weight *= cset->evaluate({this_flav == 4 ? getSystString_BTV(Correction::variation::nom) : syst_str,
                                          this_flav,
                                          fabs(jet.Eta()),
                                          jet.Pt(),
                                          this_score});
            }
        }
        return weight;
    }
    else if (method == JetTagging::JetTaggingSFMethod::comb or method == JetTagging::JetTaggingSFMethod::mujets)
    {
        auto cset = cset_btagging->at(this_taggerStr + "_" + JetTagging::GetJetTaggingSFMethodStr(method).Data());
        string light_str;
        if (Run == 2)
            light_str = this_taggerStr + "_incl";
        else if (Run == 3)
            light_str = this_taggerStr + "_light";
        auto cset_light = cset_btagging->at(light_str);
        float this_cut = GetBTaggingWP(tagger, wp);
        for (const auto &jet : jets)
        {
            int this_flav = 0;
            if (abs(jet.hadronFlavour()) == 5)
                this_flav = 5;
            else if (abs(jet.hadronFlavour()) == 4)
                this_flav = 4;
            float eff = GetBTaggingEff(jet.Eta(), jet.Pt(), this_flav, tagger, wp, syst);
            float sf;
            if (this_flav == 0)
                sf = cset_light->evaluate({syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
            else
                sf = cset->evaluate({syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
            if (jet.GetBTaggerResult(tagger) > this_cut)
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
        cout << "[Correction::GetBTaggingSF] method " << JetTagging::GetJetTaggingSFMethodStr(method) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float Correction::GetBTaggingR(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, std::string &processName, const variation syst, const TString &source) const
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string syst_str = getSystString_BTV(syst);

    if (source != "total")
        syst_str = getSystString_BTV(syst) + "_" + source;
    float weight = 1.;

    this_taggerStr += (string("_") + processName);


    auto cset = cset_btagging_R->at(this_taggerStr);

    for (const auto &jet : jets)
    {
        int this_flav = 0;
        if (abs(jet.hadronFlavour()) == 5)
            this_flav = 5;
        else if (abs(jet.hadronFlavour()) == 4)
            this_flav = 4;

        weight *= cset->evaluate({syst_str,
                                    this_flav,
                                    jet.Pt(),
                                    fabs(jet.Eta())});
    }
    return weight;
    
}

pair<float, float> Correction::GetCTaggingWP() const
{
    if ((int)global_wp >= 3)
    {
        cout << "[Correction::GetCTaggingWP] Workingpoint VeryTight and SuperTight are not available for C-Tagging" << endl;
        exit(ENODATA);
    }
    correction::Correction::Ref cset = cset_ctagging->at(global_taggerStr + "_wp_values");
    return make_pair(cset->evaluate({global_wpStr, "CvB"}), cset->evaluate({global_wpStr, "CvL"}));
}

pair<float, float> Correction::GetCTaggingWP(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp) const
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    if ((int)wp >= 3)
    {
        cout << "[Correction::GetCTaggingWP] Workingpoint VeryTight and SuperTight are not available for C-Tagging" << endl;
        exit(ENODATA);
    }
    correction::Correction::Ref cset = cset_ctagging->at(this_taggerStr + "_wp_values");
    return make_pair(cset->evaluate({this_wpStr, "CvB"}), cset->evaluate({this_wpStr, "CvL"}));
}

float Correction::GetCTaggingEff(const float eta, const float pt, const int flav, JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerWP wp, const variation syst)
{
    return 1.;
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
    correction::Correction::Ref cset = cset_btagging_eff->at(this_taggerStr);
    return cset->evaluate({
        getSystString_BTV(syst),
        this_wpStr,
    });
}

float Correction::GetCTaggingSF(const RVec<Jet> &jets, const JetTagging::JetFlavTagger tagger, const JetTagging::JetFlavTaggerWP wp, const JetTagging::JetTaggingSFMethod &method, const variation syst, const TString &source)
{
    if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet)
    {
        cerr << "[Correction::GetCTaggingSF] Run2 only supports DeepJet tagger" << endl;
        return 1.;
    }
    if (Run == 3)
    {
        //cerr << "[Correction::GetCTaggingSF] Run3 C-Tagger SF not provided by POG yet" << endl;
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
            weight *= cset->evaluate({syst_str,
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
                sf = cset->evaluate({syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
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
        cout << "[Correction::GetBTaggingSF] method " << JetTagging::GetJetTaggingSFMethodStr(method) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float Correction::GetCTaggingR(const float nTrueInt, const float HT, const JetTagging::JetFlavTagger tagger, const TString &processName, const TString &ttBarCategory, const TString &syst_str) const
{
    string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
    if (processName != "")
        this_taggerStr += "_" + processName;
    else
        this_taggerStr += (string("_") + Sample.Data());
    auto cset = cset_ctagging_R->at(this_taggerStr);
    return cset->evaluate({syst_str.Data(), ttBarCategory.Data(), nTrueInt, HT});
}

float Correction::GetJER(const float eta, const float pt, const float rho) const
{
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JER_GT.at(DataEra.Data());
    cset_string.replace(cset_string.find("######"), 6, "PtResolution");
    cset = cset_jerc->at(cset_string);
    return cset->evaluate({eta, pt, rho});
}

float Correction::GetJERSF(const float eta, const float pt, const variation syst, const TString &source) const
{
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JER_GT.at(DataEra.Data());
    cset_string.replace(cset_string.find("######"), 6, "ScaleFactor");
    cset = cset_jerc->at(cset_string);
    if (Run == 3)
    {
        return cset->evaluate({eta, pt, getSystString_JME(syst)});
    }
    else if (Run == 2)
    {
        return cset->evaluate({eta, getSystString_JME(syst)});
    }
    return 1.;
}

float Correction::GetJESUncertainty(const float eta, const float pt, const variation syst, const TString &source) const
{
    int int_syst = 0;
    if (syst == variation::up)
        int_syst = 1;
    else if (syst == variation::down)
        int_syst = -1;
    else
        int_syst = 0;

    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_JES_GT.at(DataEra.Data());
    cset_string.replace(cset_string.find("######"), 6, source);
    cset = cset_jerc->at(cset_string);
    float this_factor = 1.;
    this_factor += (int_syst * cset->evaluate({eta, pt}));
    return this_factor;
}

bool Correction::IsJetVetoZone(const float eta, const float phi, TString mapCategory) const
{
    correction::Correction::Ref cset = nullptr;
    string cset_string = JME_vetomap_keys.at(DataEra.Data());
    cset = cset_jetvetomap->at(cset_string);
    if (cset->evaluate({mapCategory.Data(), eta, phi}) > 0)
        return true;

    return false;
}

void Correction::METXYCorrection(Particle &Met, const int RunNumber, const int npvs, const XYCorrection_MetType MetType)
{
    if (Run == 3)
    {
        cout << "[Correction::METXYCorrection] Run3 does not support METXYCorrection" << endl;
        return;
    }
    correction::Correction::Ref cset_pt = nullptr;
    correction::Correction::Ref cset_phi = nullptr;
    switch (MetType)
    {
    case XYCorrection_MetType::Type1PFMET:
        if (IsDATA)
        {
            cset_pt = cset_met->at("pt_metphicorr_pfmet_data");
            cset_phi = cset_met->at("phi_metphicorr_pfmet_data");
        }
        else
        {
            cset_pt = cset_met->at("pt_metphicorr_pfmet_mc");
            cset_phi = cset_met->at("phi_metphicorr_pfmet_mc");
        }
        break;
    case XYCorrection_MetType::Type1PuppiMET:
        if (IsDATA)
        {
            cset_pt = cset_met->at("pt_metphicorr_puppimet_data");
            cset_phi = cset_met->at("phi_metphicorr_puppimet_data");
        }
        else
        {
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

float Correction::GetTopPtReweight(const RVec<Gen> &gens) const
{
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
