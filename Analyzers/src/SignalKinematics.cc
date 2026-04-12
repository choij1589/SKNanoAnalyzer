#include "SignalKinematics.h"

SignalKinematics::SignalKinematics() : channel(Channel::NONE) {}

SignalKinematics::~SignalKinematics() {}

void SignalKinematics::initializeAnalyzer() {
    TriLeptonBase::initializeAnalyzer();

    if (!Run3Mu && !Run1E2Mu) {
        throw std::runtime_error("SignalKinematics requires Run3Mu or Run1E2Mu flag");
    }
}

void SignalKinematics::executeEvent() {
    Event ev = GetEvent();

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoJetVeto || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    RVec<Gen> truth = !IsDATA ? GetAllGens() : RVec<Gen>();

    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets);

    Channel selectedChannel = selectEvent(ev, recoObjects);

    if (selectedChannel != Channel::NONE) {
        WeightInfo weights = getWeights(ev, recoObjects);
        fillObjects(selectedChannel, recoObjects, weights, truth);
    }
}

SignalKinematics::RecoObjects SignalKinematics::defineObjects(Event& ev,
                                                              const RVec<Muon>& rawMuons,
                                                              const RVec<Electron>& rawElectrons,
                                                              const RVec<Jet>& rawJets) {
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Get MET and apply Type-I correction
    Particle METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select muons (veto + tight)
    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs->GetID("tight"), 10., 2.4);

    // Select electrons (veto + tight)
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("tight"), 10., 2.5);

    // Select jets
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> jets = SelectJets(allJets, "tight", 20., max_jeteta);

    // Apply Run 2 specific filters
    if (Run == 2) {
        RVec<Jet> jets_puId;
        for (const auto& j : jets) {
            if (j.PassID("loosePuId")) jets_puId.emplace_back(j);
        }
        jets = jets_puId;

        if (!RunNoJetVeto) {
            RVec<Jet> jets_vetoMap;
            for (const auto& j : jets) {
                if (PassVetoMap(j, allMuons, "jetvetomap")) jets_vetoMap.emplace_back(j);
            }
            jets = jets_vetoMap;
        }
    }

    // Veto leptons inside jets (use veto leptons for overlap removal)
    RVec<Jet> jets_vetoLep = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);

    // Select b-jets
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& j : jets_vetoLep) {
        if (j.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet) > wp) {
            bjets.emplace_back(j);
        }
    }

    RecoObjects objects;
    objects.vetoMuons = vetoMuons;
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.tightElectrons = tightElectrons;
    objects.jets = jets_vetoLep;
    objects.bjets = bjets;
    objects.METv = METv;

    return objects;
}

SignalKinematics::Channel SignalKinematics::selectEvent(Event& ev, const RecoObjects& recoObjects) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    if (Run1E2Mu) {
        // SR1E2MU: exactly 2 tight muons + 1 tight electron, no extra leptons
        if (!(tightElectrons.size() == 1 && vetoElectrons.size() == 1 &&
              tightMuons.size() == 2 && vetoMuons.size() == 2)) return Channel::NONE;

        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];

        // OS muon pair
        if (mu1.Charge() + mu2.Charge() != 0) return Channel::NONE;

        // Dimuon mass > 12 GeV
        Particle pair = mu1 + mu2;
        if (pair.M() <= 12.) return Channel::NONE;

        // At least 2 jets
        if (jets.size() < 2) return Channel::NONE;

        // At least 1 b-jet
        if (bjets.size() == 0) return Channel::NONE;

        return Channel::SR1E2MU;
    }

    if (Run3Mu) {
        // SR3MU: exactly 3 tight muons, no electrons, no extra leptons
        if (!(tightMuons.size() == 3 && vetoMuons.size() == 3 &&
              tightElectrons.size() == 0 && vetoElectrons.size() == 0)) return Channel::NONE;

        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        const Muon& mu3 = tightMuons[2];

        // Total charge = ±1
        if (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) != 1) return Channel::NONE;

        // Configure charges: mu_ss1, mu_ss2, mu_os
        Muon mu_ss1, mu_ss2, mu_os;
        if (mu1.Charge() == mu2.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu2; mu_os = mu3;
        } else if (mu1.Charge() == mu3.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu3; mu_os = mu2;
        } else {
            mu_ss1 = mu2; mu_ss2 = mu3; mu_os = mu1;
        }

        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;
        if (pair1.M() <= 12.) return Channel::NONE;
        if (pair2.M() <= 12.) return Channel::NONE;

        // At least 2 jets
        if (jets.size() < 2) return Channel::NONE;

        // At least 1 b-jet
        if (bjets.size() == 0) return Channel::NONE;

        return Channel::SR3MU;
    }

    return Channel::NONE;
}

SignalKinematics::WeightInfo SignalKinematics::getWeights(Event& ev, const RecoObjects& recoObjects) {
    WeightInfo weights;

    if (IsDATA) {
        weights.genWeight = 1.0;
        weights.prefireWeight = 1.0;
        weights.pileupWeight = 1.0;
        weights.totWeight = 1.0;
        return weights;
    }

    weights.genWeight = MCweight() * ev.GetTriggerLumi("Full");
    weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);
    weights.pileupWeight = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    weights.totWeight = weights.genWeight * weights.prefireWeight * weights.pileupWeight;

    return weights;
}

void SignalKinematics::fillObjects(Channel ch, const RecoObjects& recoObjects, const WeightInfo& weights, const RVec<Gen>& truth) {
    TString channelStr = channelToString(ch) + "/Central";
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle& METv = recoObjects.METv;
    float w = weights.totWeight;

    // Fill weights
    FillHist(channelStr + "/weights/genWeight", weights.genWeight, 1., 200, -10000, 10000.);
    FillHist(channelStr + "/weights/prefireWeight", weights.prefireWeight, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/pileupWeight", weights.pileupWeight, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/totWeight", weights.totWeight, 1., 100, -5., 5.);

    // Fill individual muon kinematics
    for (size_t i = 0; i < muons.size(); i++) {
        TString idx = TString::Format("%zu", i+1);
        const Muon& mu = muons[i];
        FillHist(channelStr + "/muons/" + idx + "/pt", mu.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/muons/" + idx + "/eta", mu.Eta(), w, 48, -2.4, 2.4);
        FillHist(channelStr + "/muons/" + idx + "/phi", mu.Phi(), w, 64, -3.2, 3.2);
        FillHist(channelStr + "/muons/" + idx + "/mass", mu.M(), w, 10, 0., 1.);
    }

    // Fill electron kinematics (SR1E2MU)
    if (ch == Channel::SR1E2MU) {
        for (size_t i = 0; i < electrons.size(); i++) {
            TString idx = TString::Format("%zu", i+1);
            const Electron& el = electrons[i];
            FillHist(channelStr + "/electrons/" + idx + "/pt", el.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/electrons/" + idx + "/eta", el.Eta(), w, 50, -2.5, 2.5);
            FillHist(channelStr + "/electrons/" + idx + "/phi", el.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/electrons/" + idx + "/mass", el.M(), w, 100, 0., 1.);
        }
    }

    // Fill jet kinematics
    FillHist(channelStr + "/jets/size", jets.size(), w, 20, 0., 20.);
    FillHist(channelStr + "/bjets/size", bjets.size(), w, 15, 0., 15.);
    for (size_t i = 0; i < std::min(jets.size(), size_t(4)); i++) {
        TString idx = TString::Format("%zu", i+1);
        const Jet& j = jets[i];
        FillHist(channelStr + "/jets/" + idx + "/pt", j.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/jets/" + idx + "/eta", j.Eta(), w, 50, -2.5, 2.5);
    }

    // Fill MET
    FillHist(channelStr + "/METv/pt", METv.Pt(), w, 300, 0., 300.);

    // Channel-specific pair filling
    if (ch == Channel::SR1E2MU) {
        // Unambiguous OS muon pair (A -> mumu candidate)
        const Muon& mu1 = muons[0];
        const Muon& mu2 = muons[1];
        Particle dimuon = mu1 + mu2;
        FillHist(channelStr + "/dimuon/mass", dimuon.M(), w, 200, 0., 200.);
        FillHist(channelStr + "/dimuon/pt", dimuon.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/dimuon/eta", dimuon.Eta(), w, 100, -5., 5.);
        FillHist(channelStr + "/dimuon/deltaR", mu1.DeltaR(mu2), w, 100, 0., 5.);
    }

    if (ch == Channel::SR3MU) {
        // Configure charge to get two OS pairs
        const Muon& mu1 = muons[0];
        const Muon& mu2 = muons[1];
        const Muon& mu3 = muons[2];

        Muon mu_ss1, mu_ss2, mu_os;
        if (mu1.Charge() == mu2.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu2; mu_os = mu3;
        } else if (mu1.Charge() == mu3.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu3; mu_os = mu2;
        } else {
            mu_ss1 = mu2; mu_ss2 = mu3; mu_os = mu1;
        }

        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;

        FillHist(channelStr + "/pair1/mass", pair1.M(), w, 200, 0., 200.);
        FillHist(channelStr + "/pair1/pt", pair1.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/pair1/deltaR", mu_ss1.DeltaR(mu_os), w, 100, 0., 5.);
        FillHist(channelStr + "/pair2/mass", pair2.M(), w, 200, 0., 200.);
        FillHist(channelStr + "/pair2/pt", pair2.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/pair2/deltaR", mu_ss2.DeltaR(mu_os), w, 100, 0., 5.);
    }

    // Pair selection study (SR3MU, MC only)
    if (ch == Channel::SR3MU && !IsDATA && truth.size() > 0) {
        fillPairSelectionStudy(channelStr, recoObjects, truth, w);
    }

    // Gen-level distributions (MC only)
    if (!IsDATA && truth.size() > 0) {
        fillGenLeptonOrigin(channelStr, ch, recoObjects, truth, w);
        fillGenBJetMatching(channelStr, recoObjects, truth, w);
        fillGenDistributions(truth, w);
    }
}

void SignalKinematics::fillGenLeptonOrigin(const TString& channelStr, Channel ch,
                                            const RecoObjects& recoObjects,
                                            const RVec<Gen>& truth, float w) {
    if (ch == Channel::SR1E2MU) {
        const Muon& mu1 = recoObjects.tightMuons[0];
        const Muon& mu2 = recoObjects.tightMuons[1];
        const Electron& el = recoObjects.tightElectrons[0];

        // Dimuon mass from A->mumu (both muons should be type 2 = BSM-Prompt)
        int type_mu1 = GetLeptonType(mu1, truth);
        int type_mu2 = GetLeptonType(mu2, truth);
        Particle dimuon = mu1 + mu2;
        if (type_mu1 == 2 && type_mu2 == 2) {
            FillHist(channelStr + "/GenMatched/dimuon_mass_fromA", dimuon.M(), w, 200, 0., 200.);
            // Leading/subleading muons from A (already pT-ordered)
            FillHist(channelStr + "/GenMatched/muons_fromA/1/pt", mu1.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/eta", mu1.Eta(), w, 48, -2.4, 2.4);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/phi", mu1.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/mass", mu1.M(), w, 10, 0., 1.);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/pt", mu2.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/eta", mu2.Eta(), w, 48, -2.4, 2.4);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/phi", mu2.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/mass", mu2.M(), w, 10, 0., 1.);
        }
        FillHist(channelStr + "/GenMatched/muon1_type", type_mu1, w, 15, -7., 8.);
        FillHist(channelStr + "/GenMatched/muon2_type", type_mu2, w, 15, -7., 8.);

        // Electron by origin
        int type_el = GetLeptonType(el, truth);
        FillHist(channelStr + "/GenMatched/electron_type", type_el, w, 15, -7., 8.);
        if (type_el == 1) {
            // From t -> Wb -> (enu)b
            FillHist(channelStr + "/GenMatched/electron_fromW/pt", el.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/electron_fromW/eta", el.Eta(), w, 50, -2.5, 2.5);
            FillHist(channelStr + "/GenMatched/electron_fromW/phi", el.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/electron_fromW/mass", el.M(), w, 100, 0., 1.);
        } else if (type_el == 6) {
            // From t -> H+b -> WAb -> (enu)Ab, W can be offshell
            FillHist(channelStr + "/GenMatched/electron_fromOffshellW/pt", el.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/electron_fromOffshellW/eta", el.Eta(), w, 50, -2.5, 2.5);
            FillHist(channelStr + "/GenMatched/electron_fromOffshellW/phi", el.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/electron_fromOffshellW/mass", el.M(), w, 100, 0., 1.);
        }
    }

    if (ch == Channel::SR3MU) {
        const Muon& mu1 = recoObjects.tightMuons[0];
        const Muon& mu2 = recoObjects.tightMuons[1];
        const Muon& mu3 = recoObjects.tightMuons[2];

        // Get types for all muons
        int type1 = GetLeptonType(mu1, truth);
        int type2 = GetLeptonType(mu2, truth);
        int type3 = GetLeptonType(mu3, truth);

        FillHist(channelStr + "/GenMatched/muon1_type", type1, w, 15, -7., 8.);
        FillHist(channelStr + "/GenMatched/muon2_type", type2, w, 15, -7., 8.);
        FillHist(channelStr + "/GenMatched/muon3_type", type3, w, 15, -7., 8.);

        // Configure charges
        Muon mu_ss1, mu_ss2, mu_os;
        int type_ss1, type_ss2, type_os;
        if (mu1.Charge() == mu2.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu2; mu_os = mu3;
            type_ss1 = type1; type_ss2 = type2; type_os = type3;
        } else if (mu1.Charge() == mu3.Charge()) {
            mu_ss1 = mu1; mu_ss2 = mu3; mu_os = mu2;
            type_ss1 = type1; type_ss2 = type3; type_os = type2;
        } else {
            mu_ss1 = mu2; mu_ss2 = mu3; mu_os = mu1;
            type_ss1 = type2; type_ss2 = type3; type_os = type1;
        }

        // Identify signal pair (both muons type 2 from A->mumu)
        bool isSignalPair1 = (type_ss1 == 2 && type_os == 2);
        bool isSignalPair2 = (type_ss2 == 2 && type_os == 2);

        // Pair selection correctness
        // Convention: pick pair with smaller mass as signal candidate
        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;
        bool selectedPair1 = (pair1.M() < pair2.M());
        if (isSignalPair1 != isSignalPair2) {
            bool correct = (selectedPair1 == isSignalPair1);
            FillHist(channelStr + "/GenMatched/pair_selection_correct", correct ? 1.0 : 0.0, w, 2, 0., 2.);
        }

        // Fill OS dimuon mass from A->mumu and leading/subleading muons from A
        auto fillSignalPair = [&](const Muon& muA1, const Muon& muA2) {
            Particle pairA = muA1 + muA2;
            FillHist(channelStr + "/GenMatched/dimuon_mass_fromA", pairA.M(), w, 200, 0., 200.);
            const Muon& lead = (muA1.Pt() > muA2.Pt()) ? muA1 : muA2;
            const Muon& sub  = (muA1.Pt() > muA2.Pt()) ? muA2 : muA1;
            FillHist(channelStr + "/GenMatched/muons_fromA/1/pt", lead.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/eta", lead.Eta(), w, 48, -2.4, 2.4);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/phi", lead.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/muons_fromA/1/mass", lead.M(), w, 10, 0., 1.);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/pt", sub.Pt(), w, 300, 0., 300.);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/eta", sub.Eta(), w, 48, -2.4, 2.4);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/phi", sub.Phi(), w, 64, -3.2, 3.2);
            FillHist(channelStr + "/GenMatched/muons_fromA/2/mass", sub.M(), w, 10, 0., 1.);
        };
        if (isSignalPair1) fillSignalPair(mu_ss1, mu_os);
        if (isSignalPair2) fillSignalPair(mu_ss2, mu_os);

        // Muon not from A: classify by origin
        // If pair1 is signal, the other muon is mu_ss2; if pair2 is signal, it's mu_ss1
        auto fillMuonByOrigin = [&](const Muon& mu, int type) {
            if (type == 1) {
                FillHist(channelStr + "/GenMatched/muon_fromW/pt", mu.Pt(), w, 300, 0., 300.);
                FillHist(channelStr + "/GenMatched/muon_fromW/eta", mu.Eta(), w, 48, -2.4, 2.4);
                FillHist(channelStr + "/GenMatched/muon_fromW/phi", mu.Phi(), w, 64, -3.2, 3.2);
                FillHist(channelStr + "/GenMatched/muon_fromW/mass", mu.M(), w, 10, 0., 1.);
            } else if (type == 6) {
                FillHist(channelStr + "/GenMatched/muon_fromOffshellW/pt", mu.Pt(), w, 300, 0., 300.);
                FillHist(channelStr + "/GenMatched/muon_fromOffshellW/eta", mu.Eta(), w, 48, -2.4, 2.4);
                FillHist(channelStr + "/GenMatched/muon_fromOffshellW/phi", mu.Phi(), w, 64, -3.2, 3.2);
                FillHist(channelStr + "/GenMatched/muon_fromOffshellW/mass", mu.M(), w, 10, 0., 1.);
            }
        };
        if (isSignalPair1) fillMuonByOrigin(mu_ss2, type_ss2);
        if (isSignalPair2) fillMuonByOrigin(mu_ss1, type_ss1);
    }
}

void SignalKinematics::fillGenBJetMatching(const TString& channelStr,
                                            const RecoObjects& recoObjects,
                                            const RVec<Gen>& truth, float w) {
    // Find gen b-quarks with isLastCopy() and classify by top decay sibling
    struct GenBInfo {
        int index;
        float pt, eta, phi, mass;
        enum Origin { FROM_TOP_W, FROM_TOP_HPLUS, OTHER } origin;
    };

    RVec<GenBInfo> genBQuarks;
    for (size_t i = 0; i < truth.size(); i++) {
        const Gen& gen = truth[i];
        if (abs(gen.PID()) != 5) continue;
        if (!gen.isLastCopy()) continue;

        // Walk up mother chain to find the top ancestor (skip FSR b-quark copies)
        int topIdx = -1;
        int motherIdx = gen.MotherIndex();
        while (motherIdx >= 0 && motherIdx < (int)truth.size()) {
            if (abs(truth[motherIdx].PID()) == 6) {
                topIdx = motherIdx;
                break;
            }
            motherIdx = truth[motherIdx].MotherIndex();
        }

        GenBInfo::Origin origin = GenBInfo::OTHER;
        if (topIdx >= 0) {
            // Check siblings of this b-quark from the same top to distinguish t->Wb vs t->H+b
            bool hasWSibling = false;
            bool hasHplusSibling = false;
            for (size_t j = 0; j < truth.size(); j++) {
                if (j == i) continue;
                if (truth[j].MotherIndex() == topIdx) {
                    if (abs(truth[j].PID()) == 24) hasWSibling = true;
                    if (abs(truth[j].PID()) == 37) hasHplusSibling = true;
                }
            }
            if (hasHplusSibling) origin = GenBInfo::FROM_TOP_HPLUS;
            else if (hasWSibling) origin = GenBInfo::FROM_TOP_W;
        }

        GenBInfo bq;
        bq.index = i;
        bq.pt = gen.Pt();
        bq.eta = gen.Eta();
        bq.phi = gen.Phi();
        bq.mass = gen.M();
        bq.origin = origin;
        genBQuarks.push_back(bq);
    }

    // Match reco b-jets to gen b-quarks by deltaR < 0.4
    for (const auto& bjet : recoObjects.bjets) {
        float minDR = 999.;
        int bestIdx = -1;
        for (size_t i = 0; i < genBQuarks.size(); i++) {
            float deta = bjet.Eta() - genBQuarks[i].eta;
            float dphi = TVector2::Phi_mpi_pi(bjet.Phi() - genBQuarks[i].phi);
            float dr = sqrt(deta*deta + dphi*dphi);
            if (dr < minDR) {
                minDR = dr;
                bestIdx = i;
            }
        }

        if (bestIdx >= 0 && minDR < 0.4) {
            if (genBQuarks[bestIdx].origin == GenBInfo::FROM_TOP_W) {
                FillHist(channelStr + "/GenMatched/bjet_withW/pt", bjet.Pt(), w, 300, 0., 300.);
                FillHist(channelStr + "/GenMatched/bjet_withW/eta", bjet.Eta(), w, 50, -2.5, 2.5);
                FillHist(channelStr + "/GenMatched/bjet_withW/phi", bjet.Phi(), w, 64, -3.2, 3.2);
                FillHist(channelStr + "/GenMatched/bjet_withW/mass", bjet.M(), w, 100, 0., 100.);
            } else if (genBQuarks[bestIdx].origin == GenBInfo::FROM_TOP_HPLUS) {
                FillHist(channelStr + "/GenMatched/bjet_withHplus/pt", bjet.Pt(), w, 300, 0., 300.);
                FillHist(channelStr + "/GenMatched/bjet_withHplus/eta", bjet.Eta(), w, 50, -2.5, 2.5);
                FillHist(channelStr + "/GenMatched/bjet_withHplus/phi", bjet.Phi(), w, 64, -3.2, 3.2);
                FillHist(channelStr + "/GenMatched/bjet_withHplus/mass", bjet.M(), w, 100, 0., 100.);
            }
        }
    }
}

void SignalKinematics::fillGenDistributions(const RVec<Gen>& truth, float w) {
    TString prefix = "GenLevel";

    // Gen A (pseudoscalar, PID=36)
    for (const auto& gen : truth) {
        if (abs(gen.PID()) == 36 && gen.isLastCopy()) {
            FillHist(prefix + "/A/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/A/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/A/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/A/mass", gen.M(), w, 200, 0., 200.);

            // deltaR of the two child muons from A
            RVec<const Gen*> daughters;
            int genIdx = gen.Index();
            for (const auto& d : truth) {
                if (d.MotherIndex() == genIdx && abs(d.PID()) == 13 && d.Status() == 1)
                    daughters.push_back(&d);
            }
            if (daughters.size() >= 2) {
                float deta = daughters[0]->Eta() - daughters[1]->Eta();
                float dphi = TVector2::Phi_mpi_pi(daughters[0]->Phi() - daughters[1]->Phi());
                float dr = sqrt(deta*deta + dphi*dphi);
                FillHist(prefix + "/A/deltaR", dr, w, 100, 0., 5.);
            }
        }
    }

    // Gen charged Higgs (PID=37)
    for (const auto& gen : truth) {
        if (abs(gen.PID()) == 37 && gen.isLastCopy()) {
            FillHist(prefix + "/Hplus/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/Hplus/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/Hplus/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/Hplus/mass", gen.M(), w, 300, 0., 300.);
        }
    }

    // Gen leptons classified by type
    // Collect muons from A to sort by pT for leading/subleading
    RVec<const Gen*> muonsFromA;
    for (const auto& gen : truth) {
        if (gen.Status() != 1) continue;
        int absPID = abs(gen.PID());
        if (absPID != 11 && absPID != 13) continue;

        int type = GetLeptonType(gen, truth);
        TString lepName = (absPID == 11) ? "electron" : "muon";

        if (type == 2) {
            // From A (BSM-Prompt) — muons collected for pT-sorted fill below
            if (absPID == 13) muonsFromA.push_back(&gen);
        } else if (type == 1) {
            // From W (EW-Prompt)
            FillHist(prefix + "/" + lepName + "_fromW/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/" + lepName + "_fromW/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/" + lepName + "_fromW/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/" + lepName + "_fromW/mass", gen.M(), w, 10, 0., 1.);
        } else if (type == 6) {
            // From offshell W in H+ decay
            FillHist(prefix + "/" + lepName + "_fromOffshellW/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/" + lepName + "_fromOffshellW/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/" + lepName + "_fromOffshellW/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/" + lepName + "_fromOffshellW/mass", gen.M(), w, 10, 0., 1.);
        }
    }

    // Fill leading/subleading muons from A (pT-sorted)
    if (muonsFromA.size() >= 2) {
        sort(muonsFromA.begin(), muonsFromA.end(), [](const Gen* a, const Gen* b) { return a->Pt() > b->Pt(); });
        FillHist(prefix + "/muons_fromA/1/pt", muonsFromA[0]->Pt(), w, 300, 0., 300.);
        FillHist(prefix + "/muons_fromA/1/eta", muonsFromA[0]->Eta(), w, 100, -5., 5.);
        FillHist(prefix + "/muons_fromA/1/phi", muonsFromA[0]->Phi(), w, 64, -3.2, 3.2);
        FillHist(prefix + "/muons_fromA/1/mass", muonsFromA[0]->M(), w, 10, 0., 1.);
        FillHist(prefix + "/muons_fromA/2/pt", muonsFromA[1]->Pt(), w, 300, 0., 300.);
        FillHist(prefix + "/muons_fromA/2/eta", muonsFromA[1]->Eta(), w, 100, -5., 5.);
        FillHist(prefix + "/muons_fromA/2/phi", muonsFromA[1]->Phi(), w, 64, -3.2, 3.2);
        FillHist(prefix + "/muons_fromA/2/mass", muonsFromA[1]->M(), w, 10, 0., 1.);
    }

    // Gen b-quarks classified by top decay
    for (size_t i = 0; i < truth.size(); i++) {
        const Gen& gen = truth[i];
        if (abs(gen.PID()) != 5) continue;
        if (!gen.isLastCopy()) continue;

        // Walk up mother chain to find the top ancestor
        int topIdx = -1;
        int motherIdx = gen.MotherIndex();
        while (motherIdx >= 0 && motherIdx < (int)truth.size()) {
            if (abs(truth[motherIdx].PID()) == 6) {
                topIdx = motherIdx;
                break;
            }
            motherIdx = truth[motherIdx].MotherIndex();
        }
        if (topIdx < 0) continue;

        // Check siblings from the same top
        bool hasWSibling = false;
        bool hasHplusSibling = false;
        for (size_t j = 0; j < truth.size(); j++) {
            if (j == i) continue;
            if (truth[j].MotherIndex() == topIdx) {
                if (abs(truth[j].PID()) == 24) hasWSibling = true;
                if (abs(truth[j].PID()) == 37) hasHplusSibling = true;
            }
        }

        if (hasHplusSibling) {
            FillHist(prefix + "/bquark_withHplus/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/bquark_withHplus/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/bquark_withHplus/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/bquark_withHplus/mass", gen.M(), w, 100, 0., 100.);
        } else if (hasWSibling) {
            FillHist(prefix + "/bquark_withW/pt", gen.Pt(), w, 300, 0., 300.);
            FillHist(prefix + "/bquark_withW/eta", gen.Eta(), w, 100, -5., 5.);
            FillHist(prefix + "/bquark_withW/phi", gen.Phi(), w, 64, -3.2, 3.2);
            FillHist(prefix + "/bquark_withW/mass", gen.M(), w, 100, 0., 100.);
        }
    }
}

void SignalKinematics::fillPairSelectionStudy(const TString& channelStr,
                                               const RecoObjects& recoObjects,
                                               const RVec<Gen>& truth, float w) {
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const Particle& METv = recoObjects.METv;
    TString prefix = channelStr + "/PairStudy";

    // Configure SS/OS charges
    const Muon& mu1 = muons[0];
    const Muon& mu2 = muons[1];
    const Muon& mu3 = muons[2];

    Muon mu_ss1, mu_ss2, mu_os;
    int type_ss1, type_ss2, type_os;
    int type1 = GetLeptonType(mu1, truth);
    int type2 = GetLeptonType(mu2, truth);
    int type3 = GetLeptonType(mu3, truth);

    if (mu1.Charge() == mu2.Charge()) {
        mu_ss1 = mu1; mu_ss2 = mu2; mu_os = mu3;
        type_ss1 = type1; type_ss2 = type2; type_os = type3;
    } else if (mu1.Charge() == mu3.Charge()) {
        mu_ss1 = mu1; mu_ss2 = mu3; mu_os = mu2;
        type_ss1 = type1; type_ss2 = type3; type_os = type2;
    } else {
        mu_ss1 = mu2; mu_ss2 = mu3; mu_os = mu1;
        type_ss1 = type2; type_ss2 = type3; type_os = type1;
    }

    // Build pairs
    Particle pair1 = mu_ss1 + mu_os;
    Particle pair2 = mu_ss2 + mu_os;

    // Compute variables for each pair
    float mass1 = pair1.M();
    float mass2 = pair2.M();
    float gamma1 = pair1.Pt() / pair1.M();
    float gamma2 = pair2.Pt() / pair2.M();

    // mT: transverse mass of MET and the muon NOT in the pair
    // pair1 = mu_ss1 + mu_os, so the other muon is mu_ss2
    // pair2 = mu_ss2 + mu_os, so the other muon is mu_ss1
    float dphi1 = TVector2::Phi_mpi_pi(mu_ss2.Phi() - METv.Phi());
    float mT1 = sqrt(2. * mu_ss2.Pt() * METv.Pt() * (1. - cos(dphi1)));
    float dphi2 = TVector2::Phi_mpi_pi(mu_ss1.Phi() - METv.Phi());
    float mT2 = sqrt(2. * mu_ss1.Pt() * METv.Pt() * (1. - cos(dphi2)));

    // Gen truth: which pair is the signal (A->mumu)?
    bool isSignalPair1 = (type_ss1 == 2 && type_os == 2);
    bool isSignalPair2 = (type_ss2 == 2 && type_os == 2);

    // Only fill if exactly one pair is the signal (ambiguous cases skipped)
    if (!(isSignalPair1 != isSignalPair2)) return;

    // Fill signal and wrong pair distributions
    float signalMass  = isSignalPair1 ? mass1 : mass2;
    float wrongMass   = isSignalPair1 ? mass2 : mass1;
    float signalGamma = isSignalPair1 ? gamma1 : gamma2;
    float wrongGamma  = isSignalPair1 ? gamma2 : gamma1;
    float signalMT    = isSignalPair1 ? mT1 : mT2;
    float wrongMT     = isSignalPair1 ? mT2 : mT1;

    FillHist(prefix + "/signal_pair/mass",  signalMass,  w, 200, 0., 200.);
    FillHist(prefix + "/signal_pair/gamma", signalGamma, w, 100, 0., 10.);
    FillHist(prefix + "/signal_pair/mT",    signalMT,    w, 200, 0., 200.);
    FillHist(prefix + "/wrong_pair/mass",   wrongMass,   w, 200, 0., 200.);
    FillHist(prefix + "/wrong_pair/gamma",  wrongGamma,  w, 100, 0., 10.);
    FillHist(prefix + "/wrong_pair/mT",     wrongMT,     w, 200, 0., 200.);

    // Correctness: does picking the smaller value select the signal pair?
    bool massCorrect  = (signalMass < wrongMass);
    bool gammaCorrect = (signalGamma < wrongGamma);
    bool mTCorrect    = (signalMT < wrongMT);

    FillHist(prefix + "/correctness/mass",  massCorrect  ? 1.0 : 0.0, w, 2, 0., 2.);
    FillHist(prefix + "/correctness/gamma", gammaCorrect ? 1.0 : 0.0, w, 2, 0., 2.);
    FillHist(prefix + "/correctness/mT",    mTCorrect    ? 1.0 : 0.0, w, 2, 0., 2.);
}
