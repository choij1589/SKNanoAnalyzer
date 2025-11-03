#include "MatrixTreeProducer.h"

MatrixTreeProducer::MatrixTreeProducer() : channel(Channel::NONE) {}

MatrixTreeProducer::~MatrixTreeProducer() {}

void MatrixTreeProducer::WriteHist() {
    GetOutfile()->mkdir("tree");
    GetOutfile()->cd("tree");
    newtree->Write();
    GetOutfile()->cd();
}

void MatrixTreeProducer::initializeAnalyzer() {
    TriLeptonBase::initializeAnalyzer();

    // Load ParticleNet models
    // loadGraphNetModels();  // Disabled: C++ libtorch removed, use Python ParticleNet (PyAnalyzers/PromptTreeProducer.py)

    // Determine channel
    if (Run1E2Mu) channel = Channel::SR1E2Mu;
    if (Run3Mu) channel = Channel::SR3Mu;
    if (Run1E2Mu && Run3Mu) {
        throw std::runtime_error("Run1E2Mu and Run3Mu cannot be set at the same time");
    }
    if (channel == Channel::NONE) {
        throw std::runtime_error("Run1E2Mu or Run3Mu must be set");
    }

    // Initialize output tree
    GetOutfile()->cd();
    newtree = new TTree("Events", "Events");

    newtree->Branch("run", &RunNumber);
    newtree->Branch("event", &EventNumber);
    newtree->Branch("lumi", &LumiBlock);

    newtree->Branch("mass1", &mass1);
    newtree->Branch("mass2", &mass2);
    newtree->Branch("MT1", &MT1);
    newtree->Branch("MT2", &MT2);

    // Multi-class GraphNet score branches (dynamic creation per mass point and class)
    const std::vector<TString> massPoints = {"MHc160_MA85", "MHc130_MA90", "MHc100_MA95"};
    const std::vector<TString> classNames = {"signal", "nonprompt", "diboson", "ttZ"};

    for (const auto& massPoint : massPoints) {
        for (size_t i = 0; i < classNames.size(); ++i) {
            const TString& className = classNames[i];
            TString branchName;
            if (i == 0) {
                // First class (signal) uses masspoint name only
                branchName = "score_" + massPoint;
            } else {
                // Other classes append class name
                branchName = "score_" + massPoint + "_" + className;
            }
            newtree->Branch(branchName, &ParticleNetScores[massPoint][className]);
        }
    }

    newtree->Branch("fold", &fold);
    newtree->Branch("weight", &weight);
}

void MatrixTreeProducer::executeEvent() {
    Event ev = GetEvent();

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoVetoMap || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;

    RVec<Electron> rawElectrons = GetAllElectrons();

    // Initialize tree contents
    initTreeContents();

    // Process only Central (no systematics for matrix method)
    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets);
    Channel selectedChannel = selectEvent(ev, recoObjects);

    if (selectedChannel != Channel::NONE) {
        fillTree(selectedChannel, recoObjects);
    }
}

MatrixTreeProducer::RecoObjects MatrixTreeProducer::defineObjects(Event& ev,
                                                                   const RVec<Muon>& rawMuons,
                                                                   const RVec<Electron>& rawElectrons,
                                                                   const RVec<Jet>& rawJets) {
    // Create copies of the raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Get MET and apply Type-I correction
    Particle METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select objects with loose, tight selections
    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> looseMuons = SelectMuons(vetoMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(looseMuons, MuonIDs->GetID("tight"), 10., 2.4);

    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> looseElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("loose"), 15., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(looseElectrons, ElectronIDs->GetID("tight"), 15., 2.5);

    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> jets_selected = SelectJets(allJets, "tight", 20., max_jeteta);
    RVec<Jet> jets_vetoLep = JetsVetoLeptonInside(jets_selected, vetoElectrons, vetoMuons, 0.4);

    // Apply Run 2 specific filters and b-tagging
    RVec<Jet> jets;
    RVec<Jet> bjets;
    JetTagging::JetFlavTagger tagger = JetTagging::JetFlavTagger::DeepJet;
    float wp = myCorr->GetBTaggingWP(tagger, JetTagging::JetFlavTaggerWP::Medium);

    for (const auto& j : jets_vetoLep) {
        if (Run == 2) {
            if (!j.PassID("loosePuId")) continue;
            if (!(RunNoVetoMap || PassVetoMap(j, allMuons, "jetvetomap"))) continue;
        }
        jets.emplace_back(j);

        // b-tagging
        if (j.GetBTaggerResult(tagger) > wp) {
            bjets.emplace_back(j);
        }
    }

    RecoObjects objects;
    objects.vetoMuons = vetoMuons;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.jets = jets;
    objects.bjets = bjets;
    objects.METv = METv;

    return objects;
}

MatrixTreeProducer::Channel MatrixTreeProducer::selectEvent(Event& ev,
                                                             const RecoObjects& recoObjects) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& looseMuons = recoObjects.looseMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& looseElectrons = recoObjects.looseElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    bool is1E2Mu = (looseElectrons.size() == 1 && vetoElectrons.size() == 1 &&
                    looseMuons.size() == 2 && vetoMuons.size() == 2);
    bool is3Mu = (looseMuons.size() == 3 && vetoMuons.size() == 3 &&
                  looseElectrons.size() == 0 && vetoElectrons.size() == 0);

    // Check channel
    if (channel == Channel::SR1E2Mu && !is1E2Mu) return Channel::NONE;
    if (channel == Channel::SR3Mu && !is3Mu) return Channel::NONE;

    // Require at least one non-tight lepton
    if (channel == Channel::SR1E2Mu) {
        if (tightMuons.size() == looseMuons.size() && tightElectrons.size() == looseElectrons.size()) {
            return Channel::NONE;
        }
    }
    if (channel == Channel::SR3Mu) {
        if (tightMuons.size() == looseMuons.size()) {
            return Channel::NONE;
        }
    }

    // 1E2Mu baseline
    if (channel == Channel::SR1E2Mu) {
        if (!ev.PassTrigger(EMuTriggers)) return Channel::NONE;

        const Muon& mu1 = looseMuons.at(0);
        const Muon& mu2 = looseMuons.at(1);
        const Electron& ele = looseElectrons.at(0);

        bool passLeadMu = mu1.Pt() > 25. && ele.Pt() > 15.;
        bool passLeadEle = mu1.Pt() > 10. && ele.Pt() > 25.;
        bool passSafeCut = passLeadMu || passLeadEle;
        if (!passSafeCut) return Channel::NONE;
        if (mu1.Charge() + mu2.Charge() != 0) return Channel::NONE;

        Particle pair = mu1 + mu2;
        if (pair.M() <= 12.) return Channel::NONE;
        if (jets.size() < 2) return Channel::NONE;
        if (bjets.size() < 1) return Channel::NONE;

        return Channel::SR1E2Mu;
    }

    // 3Mu baseline
    if (channel == Channel::SR3Mu) {
        if (!ev.PassTrigger(DblMuTriggers)) return Channel::NONE;

        const Muon& mu1 = looseMuons.at(0);
        const Muon& mu2 = looseMuons.at(1);
        const Muon& mu3 = looseMuons.at(2);

        if (mu1.Pt() <= 20.) return Channel::NONE;
        if (mu2.Pt() <= 10.) return Channel::NONE;
        if (mu3.Pt() <= 10.) return Channel::NONE;
        if (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) != 1) return Channel::NONE;

        auto [pair1, pair2, MT1_temp, MT2_temp] = makePairs(looseMuons, recoObjects.METv);
        if (pair1.M() <= 12.) return Channel::NONE;
        if (pair2.M() <= 12.) return Channel::NONE;
        if (jets.size() < 2) return Channel::NONE;
        if (bjets.size() < 1) return Channel::NONE;

        return Channel::SR3Mu;
    }

    return Channel::NONE;
}

void MatrixTreeProducer::fillTree(Channel channel, const RecoObjects& recoObjects) {
    const RVec<Muon>& looseMuons = recoObjects.looseMuons;
    const RVec<Electron>& looseElectrons = recoObjects.looseElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle& METv = recoObjects.METv;

    // Calculate fake weight
    weight = GetFakeWeight(looseMuons, looseElectrons, "Central");

    // Fill mass variables
    if (channel == Channel::SR1E2Mu) {
        Particle pair = makePair(looseMuons);
        mass1 = pair.M();
        mass2 = -999.;
        MT1 = -999.;
        MT2 = -999.;
    } else if (channel == Channel::SR3Mu) {
        auto [pair1, pair2, this_MT1, this_MT2] = makePairs(looseMuons, METv);
        mass1 = pair1.M();
        mass2 = pair2.M();
        MT1 = this_MT1;
        MT2 = this_MT2;
    }

    // Calculate fold using METv and nJets
    int nJets = jets.size() + bjets.size();
    fold = calculateFold(METv, nJets);

    // Evaluate GraphNet scores
    evalScore(looseMuons, looseElectrons, jets, bjets, METv);

    // Fill tree
    newtree->Fill();
}

Particle MatrixTreeProducer::makePair(const RVec<Muon>& muons) {
    if (muons.size() != 2) {
        throw std::runtime_error("makePair requires exactly 2 muons");
    }
    return muons.at(0) + muons.at(1);
}

std::tuple<Particle, Particle, float, float> MatrixTreeProducer::makePairs(const RVec<Muon>& muons, const Particle& METv) {
    if (muons.size() != 3) {
        throw std::runtime_error("makePairs requires exactly 3 muons");
    }

    const Muon& mu1 = muons.at(0);
    const Muon& mu2 = muons.at(1);
    const Muon& mu3 = muons.at(2);

    Particle pair1, pair2;
    float MT1, MT2;
    if (mu1.Charge() == mu2.Charge()) {
        pair1 = mu1 + mu3;
        pair2 = mu2 + mu3;
        MT1 =  TMath::Sqrt(2 * mu1.Pt() * METv.Pt() * (1 - TMath::Cos(mu1.DeltaPhi(METv))));
        MT2 =  TMath::Sqrt(2 * mu2.Pt() * METv.Pt() * (1 - TMath::Cos(mu2.DeltaPhi(METv))));
    } else if (mu1.Charge() == mu3.Charge()) {
        pair1 = mu1 + mu2;
        pair2 = mu3 + mu2;
        MT1 =  TMath::Sqrt(2 * mu1.Pt() * METv.Pt() * (1 - TMath::Cos(mu1.DeltaPhi(METv))));
        MT2 =  TMath::Sqrt(2 * mu3.Pt() * METv.Pt() * (1 - TMath::Cos(mu3.DeltaPhi(METv))));
    } else {  // mu2.Charge() == mu3.Charge()
        pair1 = mu1 + mu2;
        pair2 = mu1 + mu3;
        MT1 =  TMath::Sqrt(2 * mu2.Pt() * METv.Pt() * (1 - TMath::Cos(mu2.DeltaPhi(METv))));
        MT2 =  TMath::Sqrt(2 * mu3.Pt() * METv.Pt() * (1 - TMath::Cos(mu3.DeltaPhi(METv))));
    }

    return {pair1, pair2, MT1, MT2};
}

void MatrixTreeProducer::initTreeContents() {
    const std::vector<TString> massPoints = {"MHc160_MA85", "MHc130_MA90", "MHc100_MA95"};
    const std::vector<TString> classNames = {"signal", "nonprompt", "diboson", "ttZ"};

    mass1 = -999.;
    mass2 = -999.;
    MT1 = -999.;
    MT2 = -999.;

    // Initialize all ParticleNet scores
    for (const auto& massPoint : massPoints) {
        for (const auto& className : classNames) {
            ParticleNetScores[massPoint][className] = -999.;
        }
    }

    fold = -999;
    weight = -999.;
}

void MatrixTreeProducer::evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                                    const RVec<Jet>& jets, const RVec<Jet>& bjets,
                                    const Particle& METv) {
    // C++ ParticleNet removed - use Python ParticleNet (PyAnalyzers/PromptTreeProducer.py)
    // Initialize scores to -999 (dummy values)
    const std::vector<TString> massPoints = {"MHc160_MA85", "MHc130_MA90", "MHc100_MA95"};
    const std::vector<TString> classNames = {"signal", "nonprompt", "diboson", "ttZ"};

    for (const auto& massPoint : massPoints) {
        for (const auto& className : classNames) {
            ParticleNetScores[massPoint][className] = -999.;
        }
    }
}
