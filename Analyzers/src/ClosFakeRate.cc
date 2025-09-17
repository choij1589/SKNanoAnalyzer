#include "ClosFakeRate.h"

ClosFakeRate::ClosFakeRate() {}

ClosFakeRate::~ClosFakeRate() {}

void ClosFakeRate::initializeAnalyzer() {
    TriLeptonBase::initializeAnalyzer();

    // Determine channel
    if (! (Run1E2Mu || Run3Mu)) {
        throw std::runtime_error("[ClosFakeRate::initializeAnalyzer] No channel specified");
    }
    if (Run1E2Mu && Run3Mu) {
        throw std::runtime_error("[ClosFakeRate::initializeAnalyzer] Cannot run both 1E2Mu and 3Mu channels");
    }
    if (RunSyst) {
        throw std::runtime_error("[ClosFakeRate::initializeAnalyzer] No systematic study for closure test");
    }
}

void ClosFakeRate::executeEvent() {
    Event ev = GetEvent();
    
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    
    RVec<Muon> rawMuons = GetAllMuons();
    // Only works for Run3
    if (! PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, "Central");
    Channel selectedChannel = selectEvent(ev, genParts, recoObjects, "Central");
    
    if (selectedChannel == Channel::NONE) return;
    
    WeightInfo centralWeights = getWeights(selectedChannel, ev, recoObjects, genParts, "Central");
    fillObjects(selectedChannel, recoObjects, centralWeights, "Central");
}

ClosFakeRate::RecoObjects ClosFakeRate::defineObjects(const Event& ev,
                                                     const RVec<Muon>& rawMuons,
                                                     const RVec<Electron>& rawElectrons,
                                                     const RVec<Jet>& rawJets,
                                                     const TString& syst) {
    RecoObjects objects;

    // Copy raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Get MET
    Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);
    objects.METv = ApplyTypeICorrection(METv, allJets, allElectrons, allMuons);

    // sort objects by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select muons
    objects.vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    objects.looseMuons = SelectMuons(objects.vetoMuons, MuonIDs->GetID("loose"), 10., 2.4);
    objects.tightMuons = SelectMuons(objects.looseMuons, MuonIDs->GetID("tight"), 10., 2.4);

    // Select electrons
    objects.vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    objects.looseElectrons = SelectElectrons(objects.vetoElectrons, ElectronIDs->GetID("loose"), 15., 2.5);
    objects.tightElectrons = SelectElectrons(objects.looseElectrons, ElectronIDs->GetID("tight"), 15., 2.5);

    // Select jets
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    objects.jets = SelectJets(allJets, "tight", 20., max_jeteta);
    if (Run == 2) {
        objects.jets = Filter(objects.jets, [&](const Jet& j) {
            return j.PassID("loosePuId");
        });
        objects.jets = Filter(objects.jets, [&](const Jet& j) {
            return PassVetoMap(j, objects.vetoMuons, "jetvetomap");
        });
    }
    objects.jets = JetsVetoLeptonInside(objects.jets, objects.vetoElectrons, objects.vetoMuons, 0.4);

    // Select b-jets
    const float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    objects.bjets = Filter(objects.jets, [&](const Jet& j) {
        return j.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet) > wp;
    });

    return objects;
}

ClosFakeRate::Channel ClosFakeRate::selectEvent(const Event& ev,
                                                const RVec<Gen>& truth,
                                                const RecoObjects& objects,
                                                const TString& syst) {
    
    bool is3Mu = (objects.looseMuons.size() == 3 && objects.vetoMuons.size() == 3 && 
                  objects.looseElectrons.size() == 0 && objects.vetoElectrons.size() == 0);
    bool is1E2Mu = (objects.looseMuons.size() == 2 && objects.vetoMuons.size() == 2 && 
                    objects.looseElectrons.size() == 1 && objects.vetoElectrons.size() == 1);

    if (Run1E2Mu) {
        if (!is1E2Mu) return Channel::NONE;
    }
    if (Run3Mu) {
        if (!is3Mu) return Channel::NONE;
    }

    // 1E2Mu channel
    if (Run1E2Mu) {
        if (!ev.PassTrigger(EMuTriggers)) return Channel::NONE;
        
        const Muon& mu1 = objects.looseMuons[0];
        const Muon& mu2 = objects.looseMuons[1];
        const Electron& ele = objects.looseElectrons[0];
        
        bool passLeadMu = mu1.Pt() > 25. && ele.Pt() > 15.;
        bool passLeadEle = mu1.Pt() > 10. && ele.Pt() > 25.;
        bool passSafeCut = passLeadMu || passLeadEle;
        if (!passSafeCut) return Channel::NONE;
        
        if (! (mu1.Charge() + mu2.Charge() == 0)) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (! (pair.M() > 12.)) return Channel::NONE;
        if (! (objects.jets.size() >= 2)) return Channel::NONE;
        if (! (objects.bjets.size() >= 1)) return Channel::NONE;
        
        if (objects.looseMuons.size() == objects.tightMuons.size() && 
            objects.looseElectrons.size() == objects.tightElectrons.size()) {
            return Channel::SR1E2MU;
        } else {
            return Channel::SB1E2MU;
        }
    }

    // 3Mu channel
    if (Run3Mu) {
        if (!ev.PassTrigger(DblMuTriggers)) return Channel::NONE;
        
        const Muon& mu1 = objects.looseMuons[0];
        const Muon& mu2 = objects.looseMuons[1];
        const Muon& mu3 = objects.looseMuons[2];
        
        if (! (mu1.Pt() > 20.)) return Channel::NONE;
        if (! (mu2.Pt() > 10.)) return Channel::NONE;
        
        if (! (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1)) return Channel::NONE;
        
        auto [mu_ss1, mu_ss2, mu_os] = configureChargeOf(objects.looseMuons);
        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;
        if (! (pair1.M() > 12.)) return Channel::NONE;
        if (! (pair2.M() > 12.)) return Channel::NONE;
        if (! (objects.jets.size() >= 2)) return Channel::NONE;
        if (! (objects.bjets.size() >= 1)) return Channel::NONE;
        
        if (objects.looseMuons.size() == objects.tightMuons.size() && 
            objects.looseElectrons.size() == objects.tightElectrons.size()) {
            return Channel::SR3MU;
        } else {
            return Channel::SB3MU;
        }
    }

    return Channel::NONE;
}

ClosFakeRate::WeightInfo ClosFakeRate::getWeights(const Channel selectedChannel,
                                                  const Event& ev,
                                                  const RecoObjects& objects,
                                                  const RVec<Gen>& genParts,
                                                  const TString& syst) {
    WeightInfo weights;
    
    weights.genWeight = IsDATA ? 1.0 : MCweight()*ev.GetTriggerLumi("Full");
    weights.prefireWeight = IsDATA ? 1.0 : GetL1PrefireWeight(MyCorrection::variation::nom);
    weights.pileupWeight = IsDATA ? 1.0 : myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    
    // Calculate fake weight for sideband regions
    if (selectedChannel == Channel::SB1E2MU || selectedChannel == Channel::SB3MU) {
        weights.fakeWeight = GetFakeWeight(objects.looseMuons, objects.looseElectrons, "QCD");
    } else {
        weights.fakeWeight = 1.;
    }
    
    weights.totalWeight = weights.genWeight * weights.prefireWeight * weights.pileupWeight * weights.fakeWeight;
    
    return weights;
}

void ClosFakeRate::fillObjects(const Channel selectedChannel,
                               const RecoObjects& objects,
                               const WeightInfo& weights,
                               const TString& syst) {
    
    TString channelStr = channelToString(selectedChannel);
    float weight = weights.totalWeight;
    
    // Fill muon histograms
    for (size_t idx = 0; idx < objects.looseMuons.size(); ++idx) {
        const Muon& mu = objects.looseMuons[idx];
        FillHist(Form("%s/%s/muons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), mu.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/muons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/muons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/muons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), mu.M(), weight, 10, 0., 1.);
    }

    // Fill electron histograms
    for (size_t idx = 0; idx < objects.looseElectrons.size(); ++idx) {
        const Electron& ele = objects.looseElectrons[idx];
        FillHist(Form("%s/%s/electrons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), ele.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/electrons/%zu/scEta", channelStr.Data(), syst.Data(), idx+1), ele.scEta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/electrons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/electrons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), ele.M(), weight, 100, 0., 1.);
    }

    // Fill jet histograms
    for (size_t idx = 0; idx < objects.jets.size(); ++idx) {
        const Jet& jet = objects.jets[idx];
        FillHist(Form("%s/%s/jets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), jet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), jet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/jets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), jet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/jets/%zu/mass", channelStr.Data(), syst.Data(), idx+1), jet.M(), weight, 100, 0., 100.);
    }

    // Fill b-jet histograms
    for (size_t idx = 0; idx < objects.bjets.size(); ++idx) {
        const Jet& bjet = objects.bjets[idx];
        FillHist(Form("%s/%s/bjets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), bjet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/bjets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), bjet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/bjets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), bjet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/bjets/%zu/mass", channelStr.Data(), syst.Data(), idx+1), bjet.M(), weight, 100, 0., 100.);
    }

    // Fill multiplicities
    FillHist(Form("%s/%s/jets/size", channelStr.Data(), syst.Data()), objects.jets.size(), weight, 20, 0., 20.);
    FillHist(Form("%s/%s/bjets/size", channelStr.Data(), syst.Data()), objects.bjets.size(), weight, 15, 0., 15.);

    // Fill MET
    FillHist(Form("%s/%s/METv/pt", channelStr.Data(), syst.Data()), objects.METv.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv/phi", channelStr.Data(), syst.Data()), objects.METv.Phi(), weight, 64, -3.2, 3.2);

    // Channel-specific histograms
    if (selectedChannel == Channel::SR1E2MU || selectedChannel == Channel::SB1E2MU) {
        Particle pair = objects.looseMuons[0] + objects.looseMuons[1];
        const Electron& nonprompt = objects.looseElectrons[0];
        
        FillHist(Form("%s/%s/pair/pt", channelStr.Data(), syst.Data()), pair.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/pair/eta", channelStr.Data(), syst.Data()), pair.Eta(), weight, 100, -5., 5.);
        FillHist(Form("%s/%s/pair/phi", channelStr.Data(), syst.Data()), pair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/pair/mass", channelStr.Data(), syst.Data()), pair.M(), weight, 200, 0., 200.);
        FillHist(Form("%s/%s/nonprompt/pt", channelStr.Data(), syst.Data()), nonprompt.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/nonprompt/eta", channelStr.Data(), syst.Data()), nonprompt.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/nonprompt/phi", channelStr.Data(), syst.Data()), nonprompt.Phi(), weight, 64, -3.2, 3.2);
    } else if (selectedChannel == Channel::SR3MU || selectedChannel == Channel::SB3MU) {
        auto [mu_ss1, mu_ss2, mu_os] = configureChargeOf(objects.looseMuons);
        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;
        
        float mZ = 91.2;
        Particle ZCand, nZCand;
        Muon nonprompt;
        if (abs(pair1.M() - mZ) < abs(pair2.M() - mZ)) {
            ZCand = pair1;
            nZCand = pair2;
            nonprompt = mu_ss2;
        } else {
            ZCand = pair2;
            nZCand = pair1;
            nonprompt = mu_ss1;
        }
        
        // Fill both pairs as "stack"
        FillHist(Form("%s/%s/stack/pt", channelStr.Data(), syst.Data()), pair1.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/stack/eta", channelStr.Data(), syst.Data()), pair1.Eta(), weight, 100, -5., 5.);
        FillHist(Form("%s/%s/stack/phi", channelStr.Data(), syst.Data()), pair1.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/stack/mass", channelStr.Data(), syst.Data()), pair1.M(), weight, 200, 0., 200.);
        FillHist(Form("%s/%s/stack/pt", channelStr.Data(), syst.Data()), pair2.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/stack/eta", channelStr.Data(), syst.Data()), pair2.Eta(), weight, 100, -5., 5.);
        FillHist(Form("%s/%s/stack/phi", channelStr.Data(), syst.Data()), pair2.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/stack/mass", channelStr.Data(), syst.Data()), pair2.M(), weight, 200, 0., 200.);
        FillHist(Form("%s/%s/nonprompt/pt", channelStr.Data(), syst.Data()), nonprompt.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/nonprompt/eta", channelStr.Data(), syst.Data()), nonprompt.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/nonprompt/phi", channelStr.Data(), syst.Data()), nonprompt.Phi(), weight, 64, -3.2, 3.2);
    }
}

std::tuple<Muon, Muon, Muon> ClosFakeRate::configureChargeOf(const RVec<Muon>& muons) {
    if (muons.size() != 3) {
        throw std::runtime_error("[ClosFakeRate::configureChargeOf] Wrong number of muons: " + std::to_string(muons.size()));
    }
    
    const Muon& mu1 = muons[0];
    const Muon& mu2 = muons[1];
    const Muon& mu3 = muons[2];
    
    if (mu1.Charge() == mu2.Charge()) {
        return std::make_tuple(mu1, mu2, mu3);
    } else if (mu1.Charge() == mu3.Charge()) {
        return std::make_tuple(mu1, mu3, mu2);
    } else if (mu2.Charge() == mu3.Charge()) {
        return std::make_tuple(mu2, mu3, mu1);
    } else {
        throw std::runtime_error("[ClosFakeRate::configureChargeOf] No same-sign pair found");
    }
}