#include "MeasTrigEff.h"

MeasTrigEff::MeasTrigEff() {}
MeasTrigEff::~MeasTrigEff() {}

void MeasTrigEff::initializeAnalyzer() {
    // Userflags
    MeasMuLegs = HasFlag("MeasMuLegs");
    MeasElLegs = HasFlag("MeasElLegs");
    MeasDblMuPairwise = HasFlag("MeasDblMuPairwise");
    MeasEMuPairwise = HasFlag("MeasEMuPairwise");

    // ID settings
    MuonIDs = new IDContainer("HcToWATight", "HcToWALoose");
    ElectronIDs = new IDContainer("HcToWATight", ((Run==2) ? "HcToWALooseRun2": "HcToWALooseRun3"));
    if (DataEra == "2016preVFP") {
        // No mass filter for DblMuTrigs, No DZ filters for EMu
        SglMuTrigs = {"HLT_IsoMu24", 
                      "HLT_IsoTkMu24"};
        SglElTrigs = {"HLT_Ele27_WPTight_Gsf"};
        DblMuTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"};
        DblMuDZTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"};
        EMuTrigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"};
        EMuDZTrigs = {};
    } else if (DataEra == "2016postVFP") {
        // No mass filter for DblMuTrigs
        SglMuTrigs = {"HLT_IsoMu24",
                      "HLT_IsoTkMu24"};
        SglElTrigs = {"HLT_Ele27_WPTight_Gsf"};
        DblMuTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL"};
        DblMuDZTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"};
        EMuTrigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"};
        EMuDZTrigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
    } else if (DataEra == "2017") {
        SglMuTrigs = {"HLT_IsoMu27"};
        SglElTrigs = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG"};
        DblMuTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"};
        DblMuDZTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"};
        DblMuDZMTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
        EMuTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
        EMuDZTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    } else if (DataEra == "2018") {  // 2018 and Run3
        SglMuTrigs = {"HLT_IsoMu24"};
        SglElTrigs = {"HLT_Ele32_WPTight_Gsf"};
        DblMuTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"};
        DblMuDZTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"};
        DblMuDZMTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
        EMuTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
        EMuDZTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    } else { // Run3 2022 / 2023
        SglMuTrigs = {"HLT_IsoMu24"};
        SglElTrigs = {"HLT_Ele30_WPTight_Gsf"};
        DblMuTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"};
        DblMuDZTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"};
        DblMuDZMTrigs = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
        EMuTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
        EMuDZTrigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    }
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample ,IsDATA);
}

void MeasTrigEff::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    bool PassPreTrigs = false;
    if (MeasElLegs) {
        PassPreTrigs = ev.PassTrigger(SglMuTrigs);
    } else if (MeasMuLegs) {
        PassPreTrigs = ev.PassTrigger(SglElTrigs);
    } else if (MeasDblMuPairwise) {
        PassPreTrigs = ev.PassTrigger(DblMuTrigs) || ev.PassTrigger(DblMuDZTrigs) || ev.PassTrigger(DblMuDZMTrigs);
    } else if (MeasEMuPairwise) {
        PassPreTrigs = ev.PassTrigger(EMuTrigs) || ev.PassTrigger(EMuDZTrigs);
    } else {
        cerr << "No valid flag set!" << endl;
        exit(EXIT_FAILURE);
    }
    if (!PassPreTrigs) return;

    // Object Definition
    RVec<Muon> rawMuons = GetMuons("NOCUT", 5., 2.4); // Rochester correction applied
    RVec<Electron> rawElectrons = GetElectrons("NOCUT", 5., 2.5);
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);
    sort(rawElectrons.begin(), rawElectrons.end(), PtComparing);

    RVec<Muon> looseMuons = SelectMuons(rawMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(rawMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> looseElectrons = SelectElectrons(rawElectrons, ElectronIDs->GetID("loose"), 15., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(rawElectrons, ElectronIDs->GetID("tight"), 15., 2.5);
    RVec<TrigObj> trigObjs = GetAllTrigObjs();

    // Weights
    float weight = 1.;
    if (!IsDATA) {
        weight = MCweight()*ev.GetTriggerLumi("Full");
        weight *= GetL1PrefireWeight();
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
    }

    if (MeasElLegs) measEMuTrigEff_ElLegs(tightMuons, looseMuons, tightElectrons, looseElectrons, trigObjs, weight);
    if (MeasMuLegs) measEMuTrigEff_MuLegs(tightMuons, looseMuons, tightElectrons, looseElectrons, trigObjs, weight);
    if (MeasDblMuPairwise) measDblMuTrigEff_DZ(ev, tightMuons, looseMuons, tightElectrons, looseElectrons, trigObjs, weight);
    if (MeasEMuPairwise) measEMuTrigEff_DZ(ev, tightMuons, looseMuons, tightElectrons, looseElectrons, trigObjs, weight);
}

void MeasTrigEff::measEMuTrigEff_MuLegs(RVec<Muon> &tightMuons, RVec<Muon> &looseMuons, RVec<Electron> &tightElectrons, RVec<Electron> &looseElectrons, RVec<TrigObj> &trigObjs, float &weight) {
    // Event Selection
    if (! (tightMuons.size() == 1 && looseMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && looseElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    
    // Trigger object matching for electron
    if (! PassSLT(el, trigObjs)) return;
    const float offline_pt_cut = (Run == 3 ? 32: (DataYear == 2016 ? 30: 35));
    if (! (el.Pt() > offline_pt_cut)) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    // Define Binning
    RVec<float> Mu23PtBins = {5., 10., 20., 23., 25., 30., 40., 50., 100., 200.};
    RVec<float> Mu8PtBins = {5., 8., 10., 15., 20., 30., 50., 100., 200.};
    RVec<float> EtaBins = {0., 0.8, 1.6, 2.4};

    const bool passMu23Leg = PassEMT_MuLeg(mu, trigObjs, 23.f);
    const bool passMu8Leg = PassEMT_MuLeg(mu, trigObjs, 8.f);

    FillHist("TrigEff_Mu23El12_MuLeg_DENOM/Central/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu23PtBins);
    FillHist("TrigEff_Mu8El23_MuLeg_DENOM/Central/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu8PtBins);
    if (passMu23Leg) FillHist("TrigEff_Mu23El12_MuLeg_NUM/Central/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu23PtBins);
    if (passMu8Leg)  FillHist("TrigEff_Mu8El23_MuLeg_NUM/Central/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu8PtBins);

    // Syst: QCD Contamination
    if (el.Pt() > 40.) {
        FillHist("TrigEff_Mu23El12_MuLeg_DENOM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu23PtBins);
        FillHist("TrigEff_Mu8El23_MuLeg_DENOM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu8PtBins);
        if (passMu23Leg) FillHist("TrigEff_Mu23El12_MuLeg_NUM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu23PtBins);
        if (passMu8Leg)  FillHist("TrigEff_Mu8El23_MuLeg_NUM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.OriginalPt(), weight, EtaBins, Mu8PtBins);
    }
}

void MeasTrigEff::measEMuTrigEff_ElLegs(RVec<Muon> &tightMuons, RVec<Muon> &looseMuons, RVec<Electron> &tightElectrons, RVec<Electron> &looseElectrons, RVec<TrigObj> &trigObjs, float &weight) {
    // Event Selection
    if (! (tightMuons.size() == 1 && looseMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && looseElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    if (! PassSLT(mu, trigObjs)) return;
    if (! (mu.Pt() > (DataYear == 2017 ? 29: 26))) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    // Define Binning
    RVec<float> El23PtBins = {5.,10.,20.,23.,25.,30.,40.,50.,100.,200.};
    RVec<float> El12PtBins = {5.,10.,12.,15.,20.,30.,50.,100.,200.};
    RVec<float> EtaBins = {0., 0.8, 1.479, 2.5};

    const bool passEl12Leg = PassEMT_ElLeg(el, trigObjs, 12.f);
    const bool passEl23Leg = PassEMT_ElLeg(el, trigObjs, 23.f);

    FillHist("TrigEff_Mu23El12_ElLeg_DENOM/Central/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El12PtBins);
    FillHist("TrigEff_Mu8El23_ElLeg_DENOM/Central/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El23PtBins);
    if (passEl12Leg) FillHist("TrigEff_Mu23El12_ElLeg_NUM/Central/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El12PtBins);
    if (passEl23Leg) FillHist("TrigEff_Mu8El23_ElLeg_NUM/Central/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El23PtBins);

    // Syst: QCD Contamination
    if (mu.Pt() > 35. && mu.PfRelIso04() < 0.15) {
        FillHist("TrigEff_Mu23El12_ElLeg_DENOM/AltTag/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El12PtBins);
        FillHist("TrigEff_Mu8El23_ElLeg_DENOM/AltTag/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El23PtBins);
        if (passEl12Leg) FillHist("TrigEff_Mu23El12_ElLeg_NUM/AltTag/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El12PtBins);
        if (passEl23Leg) FillHist("TrigEff_Mu8El23_ElLeg_NUM/AltTag/fEta_Pt", fabs(el.scEta()), el.Pt(), weight, EtaBins, El23PtBins);
    }
}

void MeasTrigEff::measEMuTrigEff_DZ(Event &ev, RVec<Muon> &tightMuons, RVec<Muon> &looseMuons, RVec<Electron> &tightElectrons, RVec<Electron> &looseElectrons, RVec<TrigObj> &trigObjs, float &weight) {
    // No DZ filter in 2016a
    if (GetEra() == "2016preVFP") return;

    // Event Selection
    if (! (tightMuons.size() == 1 && looseMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && looseElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    if (! (mu.Pt() > 10. && el.Pt() > 15.)) return;
    if (! (mu.Pt() > 25. || el.Pt() > 25.)) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    bool fire_isoMu8El23 = false;
    bool fire_isoMu8El23DZ = false;
    bool fire_isoMu23El12 = false;
    bool fire_isoMu23El12DZ = false;

    const TString pathMu8El23 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL";
    const TString pathMu23El12 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL";
    const TString pathMu8El23DZ = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
    const TString pathMu23El12DZ = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";

    if (ev.PassTrigger(pathMu8El23) && el.Pt() > 25. && mu.Pt() >10.) {
        if (PassEMT_ElLeg(el, trigObjs, 23.f) && PassEMT_MuLeg(mu, trigObjs, 8.f)) fire_isoMu8El23 = true;
        if (fire_isoMu8El23 && ev.PassTrigger(pathMu8El23DZ)) fire_isoMu8El23DZ = true;
    }

    if (ev.PassTrigger(pathMu23El12) && el.Pt() > 15. && mu.Pt() > 25.) {
        if (PassEMT_ElLeg(el, trigObjs, 12.f) && PassEMT_MuLeg(mu, trigObjs, 23.f)) fire_isoMu23El12 = true;
        if (fire_isoMu23El12 && ev.PassTrigger(pathMu23El12DZ)) fire_isoMu23El12DZ = true;
    }

    if (fire_isoMu8El23) {
        FillHist("TrigEff_Mu8El23DZ", 0., weight, 2, 0., 2.);
        if (fire_isoMu8El23DZ) FillHist("TrigEff_Mu8El23DZ", 1., weight, 2, 0., 2.);
    }
    if (fire_isoMu23El12) {
        FillHist("TrigEff_Mu23El12DZ", 0., weight, 2, 0., 2.);
        if (fire_isoMu23El12DZ) FillHist("TrigEff_Mu23El12DZ", 1., weight, 2, 0., 2.);
    }
    if (fire_isoMu8El23 || fire_isoMu23El12) {
        FillHist("TrigEff_EMuDZ", 0., weight, 2, 0., 2.);
        if ((fire_isoMu8El23 && fire_isoMu8El23DZ) || (fire_isoMu23El12 && fire_isoMu23El12DZ)) FillHist("TrigEff_EMuDZ", 1., weight, 2, 0., 2.);
    }
}

void MeasTrigEff::measDblMuTrigEff_DZ(Event &ev, RVec<Muon> &tightMuons, RVec<Muon> &looseMuons, RVec<Electron> &tightElectrons, RVec<Electron> &looseElectrons, RVec<TrigObj> &trigObjs, float &weight) {
    // Event Selection
    if (! (tightMuons.size() == 2 && looseMuons.size() == 2)) return;
    if (! (tightElectrons.size() == 0 && looseElectrons.size() == 0)) return;
    const Muon &mu1 = tightMuons.at(0);
    const Muon &mu2 = tightMuons.at(1);
    if (! (mu1.Pt() > 20. && mu2.Pt() > 10.)) return;
    if (! (mu1.Charge()+mu2.Charge() == 0)) return;

    double Mll = (mu1+mu2).M();
    double dRll = mu1.DeltaR(mu2);

    bool fire_iso = false;
    bool fire_isoDZ = false;
    bool fire_isoDZM = false;
    bool fire_isoM = false;

    for (unsigned int i = 0; i < DblMuTrigs.size(); i++) {
        const TString &pathIso = DblMuTrigs.at(i);
        const TString &pathDZ = DblMuDZTrigs.at(i);
        const TString &pathDZM = (DataYear > 2016) ? DblMuDZMTrigs.at(i) : "";
        if (! ev.PassTrigger(pathIso)) continue;

        const bool matchIso = PassDLT(mu1, trigObjs, 17.f) && PassDLT(mu2, trigObjs, 8.f);
        const bool matchIsoDZ = matchIso && ev.PassTrigger(pathDZ);
        const bool matchIsoDZM = matchIsoDZ && ev.PassTrigger(pathDZM);
        const bool matchIsoM = matchIso && ev.PassTrigger(pathDZM);

        fire_iso = fire_iso || matchIso;
        fire_isoDZ = fire_isoDZ || matchIsoDZ;
        fire_isoDZM = fire_isoDZM || matchIsoDZM;
        fire_isoM = fire_isoM || matchIsoM;
    }

    bool fire_DZ = false;
    bool fire_DZM = false;
    for (unsigned int i = 0; i < DblMuDZTrigs.size(); i++) {
        const TString &pathDZ = DblMuDZTrigs.at(i);
        const TString &pathDZM = (DataYear > 2016) ? DblMuDZMTrigs.at(i) : "";
        if (! ev.PassTrigger(pathDZ)) continue;

        const bool matchDZ = ev.PassTrigger(pathDZ);
        const bool matchDZM = matchDZ && ev.PassTrigger(pathDZM);

        fire_DZ = fire_DZ || matchDZ;
        fire_DZM = fire_DZM || matchDZM;
    }

    bool fire_M = false;
    for (unsigned int i = 0; i < DblMuDZMTrigs.size(); i++) {
        const TString &pathDZM = DblMuDZMTrigs.at(i);
        if (! ev.PassTrigger(pathDZM)) continue;

        const bool matchDZM = ev.PassTrigger(pathDZM);
        fire_M = fire_M || matchDZM;
    }

    // Fill DZ + DZM first before the mass cut
    if (fire_DZ) {
        FillHist("TrigEff_DZ", 0., weight, 1., 0., 1.);
        FillHist("TrigEff_DZ_MllNearCut", Mll, weight, 100, 0., 10.);
        FillHist("TrigEff_DZ_Mll", Mll, weight, 19, 10., 200.);
        if (fire_DZM) {
            FillHist("TrigEff_DZM", 0., weight, 1., 0., 1.);
            FillHist("TrigEff_DZM_MllNearCut", Mll, weight, 100, 0., 10.);
            FillHist("TrigEff_DZM_Mll", Mll, weight, 19, 10., 200.);
        }
    }

    if (! (fabs(Mll-91.2) < 10.)) return;
    if (! (dRll > 0.3)) return;
    if (fire_iso) {
        FillHist("TrigEff_Iso", 0., weight, 1., 0., 1);
        FillHist("TrigEff_Iso_DZ1", fabs(mu1.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_Iso_DZ2", fabs(mu2.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_Iso_DZ", fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
    }
    if (fire_isoDZ) {
        FillHist("TrigEff_IsoDZ", 0., weight, 1., 0., 1);
        FillHist("TrigEff_IsoDZ_DZ1", fabs(mu1.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_IsoDZ_DZ2", fabs(mu2.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_IsoDZ_DZ", fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
    }
    if (fire_DZ) FillHist("TrigEff_DZ", 0., weight, 1., 0., 1.);
    if (fire_isoDZM) FillHist("TrigEff_IsoDZM", 0., weight, 1., 0., 1.);
    if (fire_isoM) {
        FillHist("TrigEff_IsoM", 0., weight, 1., 0., 1);
        FillHist("TrigEff_IsoM_DZ1", fabs(mu1.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_IsoM_DZ2", fabs(mu2.dZ()), weight, 40, 0., 0.4);
        FillHist("TrigEff_IsoM_DZ", fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
    }
    if (fire_DZM) FillHist("TrigEff_DZM", 0., weight, 1., 0., 1.);
    if (fire_M) FillHist("TrigEff_M", 0., weight, 1., 0., 1.);
}

bool MeasTrigEff::PassSLT(const Electron &el, const RVec<TrigObj> &trigObjs) {
    const float trig_pt_cut = (Run == 3 ? 30: (DataYear == 2016 ? 27 : 32));
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isElectron()) continue;
        if (! (trigObj.DeltaR(el) < 0.3)) continue;
        if (! (trigObj.hasBit(1))) continue;
        if (! (trigObj.Pt() > trig_pt_cut)) continue;
        return true;
    }
    return false;
}

bool MeasTrigEff::PassSLT(const Muon &mu, const RVec<TrigObj> &trigObjs) {
    const float trig_pt_cut = (DataYear == 2017 ? 27: 24);
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(3))) continue;
        if (! (trigObj.Pt() > trig_pt_cut)) continue;
        return true;
    }
    return false;
}

bool MeasTrigEff::PassDLT(const Muon &mu, const RVec<TrigObj> &trigObjs, float pt_cut){
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(4))) continue;
        if (! (trigObj.Pt() > pt_cut)) continue;
        return true;
    }
    return false;
}

bool MeasTrigEff::PassEMT_ElLeg(const Electron &el, const RVec<TrigObj> &trigObjs, float pt_cut){
    const int trigBit = (Run == 2) ? 5 : 6; // 5 for Run2, 6 for Run3
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isElectron()) continue;
        if (! (trigObj.DeltaR(el) < 0.3)) continue;
        if (! (trigObj.hasBit(trigBit))) continue;
        if (! (trigObj.Pt() > pt_cut)) continue;
        return true;
    }
    return false;
}

bool MeasTrigEff::PassEMT_MuLeg(const Muon &mu, const RVec<TrigObj> &trigObjs, float pt_cut){
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(5))) continue;
        if (! (trigObj.Pt() > pt_cut)) continue;
        return true;
    }
    return false;
}
