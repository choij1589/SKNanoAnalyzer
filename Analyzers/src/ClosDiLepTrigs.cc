#include "ClosDiLepTrigs.h"

ClosDiLepTrigs::ClosDiLepTrigs() : channel(Channel::NONE) {}

ClosDiLepTrigs::~ClosDiLepTrigs() {}

void ClosDiLepTrigs::initializeAnalyzer() {
    DiLeptonBase::initializeAnalyzer();

    // Determine channel - check for trilepton channels first
    if (Run1E2Mu)     channel = Channel::EMUMU;
    else if (Run3Mu)  channel = Channel::MUMUMU;
    else if (RunDiMu) channel = Channel::DIMU;
    else if (RunEMu)  channel = Channel::EMU;
}

void ClosDiLepTrigs::executeEvent() {
    Event ev = GetEvent();
    
    // Basic filters
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoVetoMap || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    
    // Define objects (central only for trigger efficiency studies)
    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets);
    
    // Event selection
    Channel selectedChannel = selectEvent(ev, recoObjects);
    if (selectedChannel == Channel::NONE) return;
    
    // Calculate weights
    float weight = getWeight(selectedChannel, ev, recoObjects);
    
    // Fill trigger efficiency histograms
    fillObjects(selectedChannel, recoObjects, weight);
}

ClosDiLepTrigs::RecoObjects ClosDiLepTrigs::defineObjects(Event& ev, 
                                                         const RVec<Muon>& rawMuons, 
                                                         const RVec<Electron>& rawElectrons, 
                                                         const RVec<Jet>& rawJets) {
    // Create copies of the raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;
    
    // Get MET
    Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);

    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Object selection
    RVec<Muon> looseMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(looseMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> looseElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(looseElectrons, ElectronIDs->GetID("tight"), 10., 2.5);
    
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", 20., max_jeteta);
    if (Run == 2) tightJets = SelectJets(tightJets, "loosePuId", 20., max_jeteta);
    RVec<Jet> tightJets_vetoLep = JetsVetoLeptonInside(tightJets, looseElectrons, looseMuons, 0.4);

    // B-jet selection
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (auto& jet : tightJets_vetoLep) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    RecoObjects objects;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.tightJets = tightJets;
    objects.tightJets_vetoLep = tightJets_vetoLep;
    objects.bjets = bjets;
    objects.METv = METv;

    return objects;
}

ClosDiLepTrigs::Channel ClosDiLepTrigs::selectEvent(Event& ev, const RecoObjects& recoObjects) {
    const RVec<Muon>& looseMuons = recoObjects.looseMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& looseElectrons = recoObjects.looseElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;

    // Channel-specific lepton requirements
    bool isDiMu = (tightMuons.size() == 2 && looseMuons.size() == 2 && 
                   tightElectrons.size() == 0 && looseElectrons.size() == 0);
    bool isEMu = (tightMuons.size() == 1 && looseMuons.size() == 1 && 
                  tightElectrons.size() == 1 && looseElectrons.size() == 1);
    bool is1E2Mu = (tightMuons.size() == 2 && looseMuons.size() == 2 &&
                    tightElectrons.size() == 1 && looseElectrons.size() == 1);
    bool is3Mu = (tightMuons.size() == 3 && looseMuons.size() == 3 &&
                  tightElectrons.size() == 0 && looseElectrons.size() == 0);

    if (channel == Channel::DIMU) {
        if (!isDiMu) return Channel::NONE;
        
        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        if (mu1.Pt() <= 20.) return Channel::NONE;
        if (mu2.Pt() <= 10.) return Channel::NONE;
        return Channel::DIMU;
    }
    else if (channel == Channel::EMU) {
        if (!isEMu) return Channel::NONE;
        
        const Muon& mu = tightMuons[0];
        const Electron& ele = tightElectrons[0];
        bool leadMu = mu.Pt() > 25. && ele.Pt() > 15.;
        bool leadEle = ele.Pt() > 25. && mu.Pt() > 10.;
        if (!(leadMu || leadEle)) return Channel::NONE;
        return Channel::EMU;
    }
    else if (channel == Channel::EMUMU) {
        if (!is1E2Mu) return Channel::NONE;
        
        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        const Electron& ele = tightElectrons[0];
        
        bool leadMu = mu1.Pt() > 25. && ele.Pt() > 15.;
        bool leadEle = mu1.Pt() > 10. && ele.Pt() > 25.;
        if (!(leadMu || leadEle)) return Channel::NONE;
        
        // OS muon pair requirement
        if (mu1.Charge() + mu2.Charge() != 0) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (pair.M() <= 12.) return Channel::NONE;
        
        // Jet requirements for signal region
        if (jets.size() < 2) return Channel::NONE;
        if (bjets.size() < 1) return Channel::NONE;
        
        return Channel::EMUMU;
    }
    else if (channel == Channel::MUMUMU) {
        if (!is3Mu) return Channel::NONE;
        
        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        const Muon& mu3 = tightMuons[2];
        
        if (mu1.Pt() <= 20.) return Channel::NONE;
        if (mu2.Pt() <= 10.) return Channel::NONE;
        if (mu3.Pt() <= 10.) return Channel::NONE;
        
        // Total charge requirement
        if (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) != 1) return Channel::NONE;
        
        // All OS pairs must have mass > 12 GeV
        Particle pair1, pair2;
        if (mu1.Charge() == mu2.Charge()) {
            pair1 = mu1 + mu3;
            pair2 = mu2 + mu3;
        } else if (mu1.Charge() == mu3.Charge()) {
            pair1 = mu1 + mu2;
            pair2 = mu2 + mu3;
        } else {
            pair1 = mu1 + mu2;
            pair2 = mu1 + mu3;
        }
        
        if (pair1.M() <= 12. || pair2.M() <= 12.) return Channel::NONE;
        
        // Jet requirements for signal region
        if (jets.size() < 2) return Channel::NONE;
        if (bjets.size() < 1) return Channel::NONE;
        
        return Channel::MUMUMU;
    }

    return Channel::NONE;
}

float ClosDiLepTrigs::getWeight(const Channel& channel, const Event& event, 
                               const RecoObjects& recoObjects) {
    float weight = 1.;
    
    if (!IsDATA) {
        Event ev = GetEvent();
        weight *= MCweight() * ev.GetTriggerLumi("Full");
        weight *= GetL1PrefireWeight(MyCorrection::variation::nom);
        weight *= myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    }

    return weight;
}

void ClosDiLepTrigs::fillObjects(const Channel& channel, const RecoObjects& recoObjects, 
                                float weight) {
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    TString channelStr = channelToString(channel);
    
    // Calculate trigger efficiencies (nominal, up, down)
    float trigWeight = 1.;
    float trigWeightUp = 1.;
    float trigWeightDown = 1.;
    
    if (channel == Channel::DIMU || channel == Channel::MUMUMU) {
        trigWeight = myCorr->GetDblMuTriggerEff(muons, IsDATA, MyCorrection::variation::nom);
        trigWeightUp = myCorr->GetDblMuTriggerEff(muons, IsDATA, MyCorrection::variation::up);
        trigWeightDown = myCorr->GetDblMuTriggerEff(muons, IsDATA, MyCorrection::variation::down);
    } else if (channel == Channel::EMU || channel == Channel::EMUMU) {
        trigWeight = myCorr->GetEMuTriggerEff(electrons, muons, IsDATA, MyCorrection::variation::nom);
        trigWeightUp = myCorr->GetEMuTriggerEff(electrons, muons, IsDATA, MyCorrection::variation::up);
        trigWeightDown = myCorr->GetEMuTriggerEff(electrons, muons, IsDATA, MyCorrection::variation::down);
    }

    // Fill trigger efficiency histograms
    FillHist("sumweight", 0., weight, 5, 0., 5.);
    FillHist("sumweight", 1., weight * trigWeight, 5, 0., 5.);
    FillHist("sumweight", 2., weight * trigWeightUp, 5, 0., 5.);
    FillHist("sumweight", 3., weight * trigWeightDown, 5, 0., 5.);
    
    // Check trigger pass for events after selection
    Event ev = GetEvent();
    if ((channel == Channel::DIMU || channel == Channel::MUMUMU) && ev.PassTrigger(DblMuTriggers)) {
        FillHist("sumweight", 4., weight, 5, 0., 5.);
    }
    if ((channel == Channel::EMU || channel == Channel::EMUMU) && ev.PassTrigger(EMuTriggers)) {
        FillHist("sumweight", 4., weight, 5, 0., 5.);
    }
}

TString ClosDiLepTrigs::channelToString(const Channel& channel) {
    switch (channel) {
        case Channel::DIMU: return "DIMU";
        case Channel::EMU: return "EMU";
        case Channel::EMUMU: return "1E2MU";
        case Channel::MUMUMU: return "3MU";
        case Channel::NONE: return "NONE";
        default: return "UNKNOWN";
    }
}
