#include "TriggerStudy.h"

TriggerStudy::TriggerStudy() {}
TriggerStudy::~TriggerStudy() {}

void TriggerStudy::initializeAnalyzer() {
    // Flags
    RunEMuTrigs = HasFlag("RunEMuTrigs");
    RunEMuTrigsWithSglMuTrigs = HasFlag("RunEMuTrigsWithSglMuTrigs");
    RunEMuTrigsWithSglElTrigs = HasFlag("RunEMuTrigsWithSglElTrigs");
    RunEMuTrigsWithDblMuTrigs = HasFlag("RunEMuTrigsWithDblMuTrigs");

    if (! (DataEra == "2018"))
        throw std::runtime_error("TriggerStudy is only available for 2018 data");

    // Lepton IDs
    MuonIDs = new IDContainer("HcToWATight", "HcToWALooseRun2");
    ElectronIDs = new IDContainer("HcToWATight", "HcToWALooseRun2");

    EMuTriggers = {
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
    };
    if (RunEMuTrigsWithSglMuTrigs)
        EMuTriggers.emplace_back("HLT_IsoMu24");
    if (RunEMuTrigsWithSglElTrigs)
        EMuTriggers.emplace_back("HLT_Ele32_WPTight_Gsf");
    if (RunEMuTrigsWithDblMuTrigs)
        EMuTriggers.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample, IsDATA);    

}

void TriggerStudy::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> allJets = GetAllJets();
    if (!PassNoiseFilter(allJets, ev)) return;

    RVec<Muon> allMuons = GetAllMuons();
    RVec<Electron> allElectrons = GetAllElectrons();

    RecoObjects recoObjects = defineObjects(ev, allMuons, allElectrons, allJets, "Central");
    GenObjects genObjects;
    genObjects.genParts = GetAllGens();
    genObjects.genJets = GetAllGenJets();
    
    Channel channel = selectEvent(ev, recoObjects, genObjects,"Central");
    if (channel == Channel::NONE) return;
    float weight = 1.;
    if (!IsDATA) {
        weight = MCweight() * ev.GetTriggerLumi("Full");
        weight *= GetL1PrefireWeight(MyCorrection::variation::nom);
        weight *= myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    }
    fillObjects(channel, recoObjects, genObjects, weight, "Central");
}

TriggerStudy::RecoObjects TriggerStudy::defineObjects(Event& ev, RVec<Muon>& allMuons, RVec<Electron>& allElectrons, RVec<Jet>& allJets, const TString& syst) {
    Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);
    METv = ApplyTypeICorrection(METv, allJets, allElectrons, allMuons);

    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("tight"), 10., 2.5);

    RVec<Jet> tightJets = SelectJets(allJets, "tight", 20., 2.5);
    tightJets = SelectJets(tightJets, "loosePuId", 20., 2.5);
    RVec<Jet> tightJets_vetoMap;
    for (const auto &jet: tightJets) {
        if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.emplace_back(jet);
    }
    tightJets = tightJets_vetoMap;
    tightJets = JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4);

    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto &jet: tightJets) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    RecoObjects recoObjects;
    recoObjects.vetoMuons = vetoMuons;
    recoObjects.tightMuons = tightMuons;
    recoObjects.vetoElectrons = vetoElectrons;
    recoObjects.tightElectrons = tightElectrons;
    recoObjects.tightJets = tightJets;
    recoObjects.bjets = bjets;
    recoObjects.METv = METv;

    return recoObjects;    
}

TriggerStudy::Channel TriggerStudy::selectEvent(Event& ev, const RecoObjects& recoObjects, const GenObjects& genObjects, const TString& syst) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    bool is1E2Mu = (tightElectrons.size() == 1 && vetoElectrons.size() == 1 && 
                    tightMuons.size() == 2 && vetoMuons.size() == 2);
    bool is3Mu = (tightMuons.size() == 3 && vetoMuons.size() == 3 &&
                  tightElectrons.size() == 0 && vetoElectrons.size() == 0);

    if (! (is1E2Mu || is3Mu)) return Channel::NONE;

    if (is1E2Mu) {
        if (! ev.PassTrigger(EMuTriggers)) return Channel::NONE;
        
        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        const Electron& el = tightElectrons[0];

        // For TTLL and DYJets sample, require at least one lepton to be nonpropmt
        if ((!IsDATA) && (MCSample.Contains("TTLL") || MCSample.Contains("DYJets"))) {
            bool isLeadMuFake = GetLeptonType(mu1, genObjects.genParts) < 0;
            bool isSubleadMuFake = GetLeptonType(mu2, genObjects.genParts) < 0;
            bool isLeadElFake = GetLeptonType(el, genObjects.genParts) < 0;
            if (! (isLeadMuFake || isSubleadMuFake || isLeadElFake)) return Channel::NONE;
        }

        if (RunEMuTrigs) {
            bool leadMu = mu1.Pt() > 25. && el.Pt() > 15.;
            bool leadEl = el.Pt() > 25. && mu1.Pt() > 10.;
            if (! (leadMu || leadEl)) return Channel::NONE;
        } else if (RunEMuTrigsWithSglMuTrigs) {
            bool leadMu = mu1.Pt() > 25. && el.Pt() > 15.;
            bool leadEl = el.Pt() > 25. && mu1.Pt() > 10.;
            bool leadSglMu = mu1.Pt() > 27.;
            if (! (leadMu || leadEl || leadSglMu)) return Channel::NONE;
        } else if (RunEMuTrigsWithSglElTrigs) {
            bool leadMu = mu1.Pt() > 25. && el.Pt() > 15.;
            bool leadEl = el.Pt() > 25. && mu1.Pt() > 10.;
            bool leadSglEl = el.Pt() > 35.;
            if (! (leadMu || leadEl || leadSglEl)) return Channel::NONE;
        } else if (RunEMuTrigsWithDblMuTrigs) {
            bool leadMu = mu1.Pt() > 25. && el.Pt() > 15.;
            bool leadEl = el.Pt() > 25. && mu1.Pt() > 10.;
            bool dblMu = mu1.Pt() > 20. && mu2.Pt() > 10.;
            if (! (leadMu || leadEl || dblMu)) return Channel::NONE;
        }

        // OS muon pair requirement
        if (! (mu1.Charge() + mu2.Charge() == 0)) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (! (pair.M() > 12.)) return Channel::NONE;

        // Jet requirements for signal region
        if (! (jets.size() >= 2)) return Channel::NONE;
        if (! (bjets.size() >= 1)) return Channel::NONE;

        return Channel::SR1E2MU;
    }
    return Channel::NONE;
}

void TriggerStudy::fillCutflow(CutStage stage, const Channel& channel, float weight, const TString& syst) {
    if (syst != "Central") return;
    TString channelStr = channelToString(channel);
    if (channelStr == "None") channelStr = "PreSel";

    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/cutflow", channelStr.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}

void TriggerStudy::fillObjects(const Channel& channel, const RecoObjects& recoObjects, const GenObjects& genObjects, float weight, const TString& syst) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle METv = recoObjects.METv;

    TString channelStr = channelToString(channel);

    // Event level histograms
    FillHist(Form("%s/%s/muons/size", channelStr.Data(), syst.Data()), vetoMuons.size(), weight, 10, 0., 10.);
    FillHist(Form("%s/%s/electrons/size", channelStr.Data(), syst.Data()), vetoElectrons.size(), weight, 10, 0., 10.);
    FillHist(Form("%s/%s/jets/size", channelStr.Data(), syst.Data()), jets.size(), weight, 10, 0., 10.);
    FillHist(Form("%s/%s/bjets/size", channelStr.Data(), syst.Data()), bjets.size(), weight, 10, 0., 10.);
    FillHist(Form("%s/%s/METv/pt", channelStr.Data(), syst.Data()), METv.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv/phi", channelStr.Data(), syst.Data()), METv.Phi(), weight, 64, -3.2, 3.2);
    
    // Fill muon histograms
    for (size_t idx = 0; idx < vetoMuons.size(); ++idx) {
        const Muon& mu = vetoMuons.at(idx);
        FillHist(Form("%s/%s/muons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), mu.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/muons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/muons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/muons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), mu.M(), weight, 10, 0., 1.);
    }
    
    // Fill electron histograms
    for (size_t idx = 0; idx < vetoElectrons.size(); ++idx) {
        const Electron& ele = vetoElectrons.at(idx);
        FillHist(Form("%s/%s/electrons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), ele.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/electrons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), ele.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/electrons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/electrons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), ele.M(), weight, 100, 0., 1.);
    }
    
    // Fill jet histograms
    for (size_t idx = 0; idx < jets.size(); ++idx) {
        const Jet& jet = jets.at(idx);
        FillHist(Form("%s/%s/jets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), jet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/rawPt", channelStr.Data(), syst.Data(), idx+1), jet.GetRawPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/originalPt", channelStr.Data(), syst.Data(), idx+1), jet.GetOriginalPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), jet.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/jets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), jet.Phi(), weight, 64, -3.2, 3.2);
    }

    // Fill bjet histograms
    for (size_t idx = 0; idx < bjets.size(); ++idx) {
        const Jet& bjet = bjets.at(idx);
        FillHist(Form("%s/%s/bjets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), bjet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/bjets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), bjet.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/bjets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), bjet.Phi(), weight, 64, -3.2, 3.2);
    }

    Particle pair = tightMuons[0] + tightMuons[1];
    FillHist(Form("%s/%s/pair/pt", channelStr.Data(), syst.Data()), pair.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/pair/eta", channelStr.Data(), syst.Data()), pair.Eta(), weight, 50, -2.5, 2.5);
    FillHist(Form("%s/%s/pair/phi", channelStr.Data(), syst.Data()), pair.Phi(), weight, 64, -3.2, 3.2);
    FillHist(Form("%s/%s/pair/mass", channelStr.Data(), syst.Data()), pair.M(), weight, 200, 0., 200.);
}
