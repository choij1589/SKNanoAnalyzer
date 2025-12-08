#include "SignalKinematics.h"

SignalKinematics::SignalKinematics() : channel(Channel::NONE) {}

SignalKinematics::~SignalKinematics() {}

void SignalKinematics::initializeAnalyzer() {
    TriLeptonBase::initializeAnalyzer();

    // This analyzer only runs on 3Mu channel
    if (!Run3Mu) {
        throw std::runtime_error("SignalKinematics requires Run3Mu flag");
    }

    channel = Channel::SR3MU;
}

void SignalKinematics::executeEvent() {
    Event ev = GetEvent();

    // Initial weight
    float initialWeight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoVetoMap || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    RVec<Gen> truth = !IsDATA ? GetAllGens() : RVec<Gen>();

    // Define objects (no systematics, simplified)
    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets);

    // Event selection
    Channel selectedChannel = selectEvent(ev, recoObjects);

    if (selectedChannel == Channel::SR3MU) {
        WeightInfo weights = getWeights(ev, recoObjects);
        fillObjects(selectedChannel, recoObjects, weights, truth);
    }
}

SignalKinematics::RecoObjects SignalKinematics::defineObjects(Event& ev,
                                                              const RVec<Muon>& rawMuons,
                                                              const RVec<Electron>& rawElectrons,
                                                              const RVec<Jet>& rawJets) {
    // Create copies
    RVec<Muon> allMuons = rawMuons;
    RVec<Jet> allJets = rawJets;

    // Get MET and apply Type-I correction
    Particle METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    Particle METv = ApplyTypeICorrection(METv_default, allJets, rawElectrons, allMuons);

    // Sort by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select tight muons
    RVec<Muon> tightMuons = SelectMuons(allMuons, MuonIDs->GetID("tight"), 10., 2.4);

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

        if (!RunNoVetoMap) {
            RVec<Jet> jets_vetoMap;
            for (const auto& j : jets) {
                if (PassVetoMap(j, allMuons, "jetvetomap")) jets_vetoMap.emplace_back(j);
            }
            jets = jets_vetoMap;
        }
    }

    // Veto leptons inside jets
    RVec<Jet> jets_vetoLep = JetsVetoLeptonInside(jets, rawElectrons, tightMuons, 0.4);

    // Select b-jets
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& j : jets_vetoLep) {
        if (j.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet) > wp) {
            bjets.emplace_back(j);
        }
    }

    RecoObjects objects;
    objects.tightMuons = tightMuons;
    objects.jets = jets_vetoLep;
    objects.bjets = bjets;
    objects.METv = METv;

    return objects;
}

SignalKinematics::Channel SignalKinematics::selectEvent(Event& ev, const RecoObjects& recoObjects) {
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    // Require exactly 3 tight muons
    if (tightMuons.size() != 3) return Channel::NONE;

    // Pass DblMu triggers
    if (!ev.PassTrigger(DblMuTriggers)) return Channel::NONE;

    const Muon& mu1 = tightMuons[0];
    const Muon& mu2 = tightMuons[1];
    const Muon& mu3 = tightMuons[2];

    // pT cuts
    if (mu1.Pt() <= 20.) return Channel::NONE;
    if (mu2.Pt() <= 10.) return Channel::NONE;
    if (mu3.Pt() <= 10.) return Channel::NONE;

    // Total charge = ±1
    if (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) != 1) return Channel::NONE;

    // Find OS pairs and check mass > 12 GeV
    // Configure charges: mu_ss1, mu_ss2, mu_os
    Muon mu_ss1, mu_ss2, mu_os;
    if (mu1.Charge() == mu2.Charge()) {
        mu_ss1 = mu1; mu_ss2 = mu2; mu_os = mu3;
    } else if (mu1.Charge() == mu3.Charge()) {
        mu_ss1 = mu1; mu_ss2 = mu3; mu_os = mu2;
    } else if (mu2.Charge() == mu3.Charge()) {
        mu_ss1 = mu2; mu_ss2 = mu3; mu_os = mu1;
    } else {
        return Channel::NONE;
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

SignalKinematics::WeightInfo SignalKinematics::getWeights(Event& ev, const RecoObjects& recoObjects) {
    WeightInfo weights;

    if (IsDATA) {
        weights.genWeight = 1.0;
        weights.prefireWeight = 1.0;
        weights.pileupWeight = 1.0;
        weights.muonRecoSF = 1.0;
        weights.muonIDSF = 1.0;
        weights.trigSF = 1.0;
        weights.btagSF = 1.0;
        weights.totWeight = 1.0;
        return weights;
    }

    weights.genWeight = MCweight() * ev.GetTriggerLumi("Full");
    weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);
    weights.pileupWeight = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    weights.muonRecoSF = myCorr->GetMuonRECOSF(recoObjects.tightMuons);
    weights.muonIDSF = myCorr->GetMuonIDSF("TopHNT", recoObjects.tightMuons);
    weights.trigSF = myCorr->GetDblMuTriggerSF(recoObjects.tightMuons);
    weights.btagSF = myCorr->GetBTaggingReweightMethod1a(
        recoObjects.jets,
        JetTagging::JetFlavTagger::DeepJet,
        JetTagging::JetFlavTaggerWP::Medium,
        JetTagging::JetTaggingSFMethod::mujets,
        MyCorrection::variation::nom,
        "central"
    );

    weights.totWeight = weights.genWeight * weights.prefireWeight * weights.pileupWeight *
                       weights.muonRecoSF * weights.muonIDSF * weights.trigSF * weights.btagSF;

    return weights;
}

void SignalKinematics::fillObjects(Channel ch, const RecoObjects& recoObjects, const WeightInfo& weights, const RVec<Gen>& truth) {
    TString channelStr = channelToString(ch);
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle& METv = recoObjects.METv;
    float w = weights.totWeight;

    // Fill weights
    FillHist(channelStr + "/weights/genWeight", weights.genWeight, 1., 200, -10000, 10000.);
    FillHist(channelStr + "/weights/prefireWeight", weights.prefireWeight, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/pileupWeight", weights.pileupWeight, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/muonRecoSF", weights.muonRecoSF, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/muonIDSF", weights.muonIDSF, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/trigSF", weights.trigSF, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/btagSF", weights.btagSF, 1., 100, -5., 5.);
    FillHist(channelStr + "/weights/totWeight", weights.totWeight, 1., 100, -5., 5.);

    // Fill individual muon kinematics
    for (size_t i = 0; i < muons.size(); i++) {
        TString idx = TString::Format("%zu", i+1);
        const Muon& mu = muons[i];
        FillHist(channelStr + "/muons/" + idx + "/pt", mu.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/muons/" + idx + "/eta", mu.Eta(), w, 48, -2.4, 2.4);
        FillHist(channelStr + "/muons/" + idx + "/phi", mu.Phi(), w, 64, -3.2, 3.2);
        FillHist(channelStr + "/muons/" + idx + "/miniIso", mu.MiniPFRelIso(), w, 100, 0., 1.);
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

    // Helper function to fill pair histograms
    auto fillPairHistograms = [&](TString prefix, const Particle& pair, const Muon& muA, const Muon& muB, const Muon& mu3rd, const RVec<Jet>& bjets, const RVec<Jet>& jets, const Particle& METv) {
        // Basic kinematics
        FillHist(channelStr + "/" + prefix + "/mass", pair.M(), w, 200, 0., 200.);
        FillHist(channelStr + "/" + prefix + "/pt", pair.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/" + prefix + "/eta", pair.Eta(), w, 100, -5., 5.);
        FillHist(channelStr + "/" + prefix + "/phi", pair.Phi(), w, 64, -3.2, 3.2);

        // Separation variables
        float deltaR = muA.DeltaR(muB);
        float deltaEta = abs(muA.Eta() - muB.Eta());
        float deltaPhi = abs(muA.DeltaPhi(muB));
        float acoplanarity = TMath::Pi() - deltaPhi;
        FillHist(channelStr + "/" + prefix + "/deltaR", deltaR, w, 100, 0., 5.);
        FillHist(channelStr + "/" + prefix + "/deltaEta", deltaEta, w, 100, 0., 5.);
        FillHist(channelStr + "/" + prefix + "/deltaPhi", deltaPhi, w, 100, 0., 5.);
        FillHist(channelStr + "/" + prefix + "/acoplanarity", acoplanarity, w, 100, 0., 3.15);

        // pT asymmetry
        float ptAsymmetry = (muA.Pt() - muB.Pt()) / (muA.Pt() + muB.Pt());
        FillHist(channelStr + "/" + prefix + "/ptAsymmetry", ptAsymmetry, w, 100, -1., 1.);

        // Scalar pT sum
        float scalarPtSum = muA.Pt() + muB.Pt();
        FillHist(channelStr + "/" + prefix + "/scalarPtSum", scalarPtSum, w, 400, 0., 400.);

        // Gamma factor (boost)
        float gammaFactor = pair.Pt() / pair.M();
        FillHist(channelStr + "/" + prefix + "/gammaFactor", gammaFactor, w, 100, 0., 10.);

        // Variables with 3rd muon
        float deltaR_pair_mu3rd = pair.DeltaR(mu3rd);
        float deltaPhi_pair_mu3rd = abs(pair.DeltaPhi(mu3rd));
        float deltaEta_pair_mu3rd = abs(pair.Eta() - mu3rd.Eta());
        float ptRatio_mu3rd = mu3rd.Pt() / pair.Pt();
        FillHist(channelStr + "/" + prefix + "/deltaR_pair_mu3rd", deltaR_pair_mu3rd, w, 100, 0., 6.);
        FillHist(channelStr + "/" + prefix + "/deltaPhi_pair_mu3rd", deltaPhi_pair_mu3rd, w, 100, 0., 3.15);
        FillHist(channelStr + "/" + prefix + "/deltaEta_pair_mu3rd", deltaEta_pair_mu3rd, w, 100, 0., 5.);
        FillHist(channelStr + "/" + prefix + "/ptRatio_mu3rd", ptRatio_mu3rd, w, 100, 0., 5.);

        // deltaR to nearest b-jet
        float deltaR_nearestBjet = 999.;
        for (const auto& bjet : bjets) {
            float dr = pair.DeltaR(bjet);
            if (dr < deltaR_nearestBjet) deltaR_nearestBjet = dr;
        }
        FillHist(channelStr + "/" + prefix + "/deltaR_nearestBjet", deltaR_nearestBjet, w, 100, 0., 6.);

        // MET-related variables
        float deltaPhi_pair_MET = abs(pair.DeltaPhi(METv));
        FillHist(channelStr + "/" + prefix + "/deltaPhi_pair_MET", deltaPhi_pair_MET, w, 100, 0., 3.15);

        float MT_pair_MET = sqrt(2 * pair.Pt() * METv.Pt() * (1 - cos(pair.DeltaPhi(METv))));
        FillHist(channelStr + "/" + prefix + "/MT_pair_MET", MT_pair_MET, w, 150, 0., 300.);

        // Same-sign muon vs MET (muA is the same-sign muon)
        float deltaPhi_muSS_MET = abs(muA.DeltaPhi(METv));
        FillHist(channelStr + "/" + prefix + "/deltaPhi_muSS_MET", deltaPhi_muSS_MET, w, 100, 0., 3.15);

        float MT_muSS_MET = sqrt(2 * muA.Pt() * METv.Pt() * (1 - cos(muA.DeltaPhi(METv))));
        FillHist(channelStr + "/" + prefix + "/MT_muSS_MET", MT_muSS_MET, w, 150, 0., 300.);

        // MT asymmetry between constituent muons
        float MT_muA_MET = sqrt(2 * muA.Pt() * METv.Pt() * (1 - cos(muA.DeltaPhi(METv))));
        float MT_muB_MET = sqrt(2 * muB.Pt() * METv.Pt() * (1 - cos(muB.DeltaPhi(METv))));
        float MT_asymmetry = abs(MT_muA_MET - MT_muB_MET) / (MT_muA_MET + MT_muB_MET);
        FillHist(channelStr + "/" + prefix + "/MT_asymmetry", MT_asymmetry, w, 100, 0., 1.);

        // deltaR to leading non-b jet
        float deltaR_leadingNonBjet = 999.;
        for (const auto& jet : jets) {
            bool isBjet = false;
            for (const auto& bjet : bjets) {
                if (jet.DeltaR(bjet) < 0.01) { isBjet = true; break; }
            }
            if (!isBjet) {
                deltaR_leadingNonBjet = pair.DeltaR(jet);
                break;  // jets are pT-sorted, so first non-b is leading
            }
        }
        FillHist(channelStr + "/" + prefix + "/deltaR_leadingNonBjet", deltaR_leadingNonBjet, w, 100, 0., 6.);

        // Individual muon properties
        FillHist(channelStr + "/" + prefix + "/mu1_pt", muA.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/" + prefix + "/mu2_pt", muB.Pt(), w, 300, 0., 300.);
        FillHist(channelStr + "/" + prefix + "/mu1_iso", muA.MiniPFRelIso(), w, 100, 0., 1.);
        FillHist(channelStr + "/" + prefix + "/mu2_iso", muB.MiniPFRelIso(), w, 100, 0., 1.);
    };

    // Truth-matching for signal vs fake discrimination
    if (!IsDATA && truth.size() > 0) {
        // Check pair1
        int type_ss1 = GetLeptonType(mu_ss1, truth);
        int type_os_for_pair1 = GetLeptonType(mu_os, truth);
        bool isSignalPair1 = (type_ss1 == 2 && type_os_for_pair1 == 2);

        if (isSignalPair1) {
            fillPairHistograms("SignalPair", pair1, mu_ss1, mu_os, mu_ss2, bjets, jets, METv);
        } else {
            fillPairHistograms("FakePair", pair1, mu_ss1, mu_os, mu_ss2, bjets, jets, METv);
        }

        // Check pair2
        int type_ss2 = GetLeptonType(mu_ss2, truth);
        int type_os_for_pair2 = GetLeptonType(mu_os, truth);
        bool isSignalPair2 = (type_ss2 == 2 && type_os_for_pair2 == 2);

        if (isSignalPair2) {
            fillPairHistograms("SignalPair", pair2, mu_ss2, mu_os, mu_ss1, bjets, jets, METv);
        } else {
            fillPairHistograms("FakePair", pair2, mu_ss2, mu_os, mu_ss1, bjets, jets, METv);
        }

        // Discrimination power analysis - only if exactly one pair is signal
        if (isSignalPair1 != isSignalPair2) {
            // Calculate discriminating variables for both pairs
            float acop1 = TMath::Pi() - abs(mu_ss1.DeltaPhi(mu_os));
            float acop2 = TMath::Pi() - abs(mu_ss2.DeltaPhi(mu_os));
            float ptSum1 = mu_ss1.Pt() + mu_os.Pt();
            float ptSum2 = mu_ss2.Pt() + mu_os.Pt();
            float ptAsym1 = abs(mu_ss1.Pt() - mu_os.Pt()) / (mu_ss1.Pt() + mu_os.Pt());
            float ptAsym2 = abs(mu_ss2.Pt() - mu_os.Pt()) / (mu_ss2.Pt() + mu_os.Pt());

            // Test acoplanarity: pick pair with smaller acoplanarity (more back-to-back)
            bool selected_pair1_acop = (acop1 < acop2);
            bool acop_correct = (selected_pair1_acop == isSignalPair1);
            FillHist(channelStr + "/Discrimination/acoplanarity_correct", acop_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test scalarPtSum: pick pair with smaller pT sum
            bool selected_pair1_ptSum = (ptSum1 < ptSum2);
            bool ptSum_correct = (selected_pair1_ptSum == isSignalPair1);
            FillHist(channelStr + "/Discrimination/scalarPtSum_correct", ptSum_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test ptAsymmetry: pick pair with smaller asymmetry (more balanced)
            bool selected_pair1_ptAsym = (ptAsym1 < ptAsym2);
            bool ptAsym_correct = (selected_pair1_ptAsym == isSignalPair1);
            FillHist(channelStr + "/Discrimination/ptAsymmetry_correct", ptAsym_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate gamma factor for both pairs
            float gamma1 = pair1.Pt() / pair1.M();
            float gamma2 = pair2.Pt() / pair2.M();

            // Test gamma factor: pick pair with smaller gamma (less boosted)
            bool selected_pair1_gamma_smaller = (gamma1 < gamma2);
            bool gamma_smaller_correct = (selected_pair1_gamma_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/gammaFactor_smaller_correct", gamma_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test gamma factor: pick pair with larger gamma (more boosted)
            bool selected_pair1_gamma_larger = (gamma1 > gamma2);
            bool gamma_larger_correct = (selected_pair1_gamma_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/gammaFactor_larger_correct", gamma_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate gamma * acoplanarity for both pairs
            float gammaAcop1 = acop1 / gamma1;
            float gammaAcop2 = acop2 / gamma2;

            // Test gamma*acop: pick pair with smaller value
            bool selected_pair1_gammaAcop_smaller = (gammaAcop1 < gammaAcop2);
            bool gammaAcop_smaller_correct = (selected_pair1_gammaAcop_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/gammaAcop_smaller_correct", gammaAcop_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test gamma*acop: pick pair with larger value
            bool selected_pair1_gammaAcop_larger = (gammaAcop1 > gammaAcop2);
            bool gammaAcop_larger_correct = (selected_pair1_gammaAcop_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/gammaAcop_larger_correct", gammaAcop_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate 3rd muon variables for both pairs
            // For pair1, 3rd muon is mu_ss2; for pair2, 3rd muon is mu_ss1
            float deltaR_pair1_mu3rd = pair1.DeltaR(mu_ss2);
            float deltaR_pair2_mu3rd = pair2.DeltaR(mu_ss1);
            float deltaPhi_pair1_mu3rd = abs(pair1.DeltaPhi(mu_ss2));
            float deltaPhi_pair2_mu3rd = abs(pair2.DeltaPhi(mu_ss1));
            float deltaEta_pair1_mu3rd = abs(pair1.Eta() - mu_ss2.Eta());
            float deltaEta_pair2_mu3rd = abs(pair2.Eta() - mu_ss1.Eta());
            float ptRatio_pair1_mu3rd = mu_ss2.Pt() / pair1.Pt();
            float ptRatio_pair2_mu3rd = mu_ss1.Pt() / pair2.Pt();

            // Test deltaR with 3rd muon: pick pair with larger separation
            bool selected_pair1_deltaR_mu3rd_larger = (deltaR_pair1_mu3rd > deltaR_pair2_mu3rd);
            bool deltaR_mu3rd_larger_correct = (selected_pair1_deltaR_mu3rd_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_pair_mu3rd_larger_correct", deltaR_mu3rd_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test deltaR with 3rd muon: pick pair with smaller separation
            bool selected_pair1_deltaR_mu3rd_smaller = (deltaR_pair1_mu3rd < deltaR_pair2_mu3rd);
            bool deltaR_mu3rd_smaller_correct = (selected_pair1_deltaR_mu3rd_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_pair_mu3rd_smaller_correct", deltaR_mu3rd_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test deltaPhi with 3rd muon: pick pair with larger phi separation
            bool selected_pair1_deltaPhi_mu3rd_larger = (deltaPhi_pair1_mu3rd > deltaPhi_pair2_mu3rd);
            bool deltaPhi_mu3rd_larger_correct = (selected_pair1_deltaPhi_mu3rd_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaPhi_pair_mu3rd_larger_correct", deltaPhi_mu3rd_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test ptRatio with 3rd muon: pick pair with smaller pT ratio
            bool selected_pair1_ptRatio_mu3rd_smaller = (ptRatio_pair1_mu3rd < ptRatio_pair2_mu3rd);
            bool ptRatio_mu3rd_smaller_correct = (selected_pair1_ptRatio_mu3rd_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/ptRatio_mu3rd_smaller_correct", ptRatio_mu3rd_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate deltaR to nearest b-jet for both pairs
            float deltaR_nearestBjet1 = 999.;
            for (const auto& bjet : bjets) {
                float dr = pair1.DeltaR(bjet);
                if (dr < deltaR_nearestBjet1) deltaR_nearestBjet1 = dr;
            }
            float deltaR_nearestBjet2 = 999.;
            for (const auto& bjet : bjets) {
                float dr = pair2.DeltaR(bjet);
                if (dr < deltaR_nearestBjet2) deltaR_nearestBjet2 = dr;
            }

            // Test deltaR to nearest b-jet: pick pair with smaller distance (closer to b-jet)
            bool selected_pair1_bjet_smaller = (deltaR_nearestBjet1 < deltaR_nearestBjet2);
            bool bjet_smaller_correct = (selected_pair1_bjet_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_nearestBjet_smaller_correct", bjet_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test deltaR to nearest b-jet: pick pair with larger distance (further from b-jet)
            bool selected_pair1_bjet_larger = (deltaR_nearestBjet1 > deltaR_nearestBjet2);
            bool bjet_larger_correct = (selected_pair1_bjet_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_nearestBjet_larger_correct", bjet_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate MET-related variables for both pairs
            float deltaPhi_pair1_MET = abs(pair1.DeltaPhi(METv));
            float deltaPhi_pair2_MET = abs(pair2.DeltaPhi(METv));
            float MT_pair1_MET = sqrt(2 * pair1.Pt() * METv.Pt() * (1 - cos(pair1.DeltaPhi(METv))));
            float MT_pair2_MET = sqrt(2 * pair2.Pt() * METv.Pt() * (1 - cos(pair2.DeltaPhi(METv))));

            // Test deltaPhi(pair, MET): pick pair with smaller deltaPhi
            bool selected_pair1_deltaPhi_MET_smaller = (deltaPhi_pair1_MET < deltaPhi_pair2_MET);
            bool deltaPhi_MET_smaller_correct = (selected_pair1_deltaPhi_MET_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaPhi_pair_MET_smaller_correct", deltaPhi_MET_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test deltaPhi(pair, MET): pick pair with larger deltaPhi
            bool selected_pair1_deltaPhi_MET_larger = (deltaPhi_pair1_MET > deltaPhi_pair2_MET);
            bool deltaPhi_MET_larger_correct = (selected_pair1_deltaPhi_MET_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaPhi_pair_MET_larger_correct", deltaPhi_MET_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test MT(pair, MET): pick pair with smaller MT
            bool selected_pair1_MT_MET_smaller = (MT_pair1_MET < MT_pair2_MET);
            bool MT_MET_smaller_correct = (selected_pair1_MT_MET_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_pair_MET_smaller_correct", MT_MET_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test MT(pair, MET): pick pair with larger MT
            bool selected_pair1_MT_MET_larger = (MT_pair1_MET > MT_pair2_MET);
            bool MT_MET_larger_correct = (selected_pair1_MT_MET_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_pair_MET_larger_correct", MT_MET_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Calculate deltaR to leading non-b jet for both pairs
            float deltaR_leadingNonBjet1 = 999.;
            for (const auto& jet : jets) {
                bool isBjet = false;
                for (const auto& bjet : bjets) {
                    if (jet.DeltaR(bjet) < 0.01) { isBjet = true; break; }
                }
                if (!isBjet) {
                    deltaR_leadingNonBjet1 = pair1.DeltaR(jet);
                    break;
                }
            }
            float deltaR_leadingNonBjet2 = 999.;
            for (const auto& jet : jets) {
                bool isBjet = false;
                for (const auto& bjet : bjets) {
                    if (jet.DeltaR(bjet) < 0.01) { isBjet = true; break; }
                }
                if (!isBjet) {
                    deltaR_leadingNonBjet2 = pair2.DeltaR(jet);
                    break;
                }
            }

            // Test deltaR to leading non-b jet: pick pair with smaller distance
            bool selected_pair1_leadingNonBjet_smaller = (deltaR_leadingNonBjet1 < deltaR_leadingNonBjet2);
            bool leadingNonBjet_smaller_correct = (selected_pair1_leadingNonBjet_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_leadingNonBjet_smaller_correct", leadingNonBjet_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test deltaR to leading non-b jet: pick pair with larger distance
            bool selected_pair1_leadingNonBjet_larger = (deltaR_leadingNonBjet1 > deltaR_leadingNonBjet2);
            bool leadingNonBjet_larger_correct = (selected_pair1_leadingNonBjet_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaR_leadingNonBjet_larger_correct", leadingNonBjet_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Same-sign muon vs MET discrimination
            // mu_ss from A (signal): LARGER deltaPhi (less aligned with MET)
            // mu_ss from W (fake): SMALLER deltaPhi (more aligned with neutrino in MET)
            float deltaPhi_ss1_MET = abs(mu_ss1.DeltaPhi(METv));
            float deltaPhi_ss2_MET = abs(mu_ss2.DeltaPhi(METv));
            float MT_ss1_MET = sqrt(2 * mu_ss1.Pt() * METv.Pt() * (1 - cos(mu_ss1.DeltaPhi(METv))));
            float MT_ss2_MET = sqrt(2 * mu_ss2.Pt() * METv.Pt() * (1 - cos(mu_ss2.DeltaPhi(METv))));

            // Test: pick mu_ss with LARGER deltaPhi (signal should be less aligned with MET)
            bool selected_pair1_deltaPhi_muSS_larger = (deltaPhi_ss1_MET > deltaPhi_ss2_MET);
            bool deltaPhi_muSS_larger_correct = (selected_pair1_deltaPhi_muSS_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaPhi_muSS_MET_larger_correct", deltaPhi_muSS_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test: pick mu_ss with SMALLER deltaPhi
            bool selected_pair1_deltaPhi_muSS_smaller = (deltaPhi_ss1_MET < deltaPhi_ss2_MET);
            bool deltaPhi_muSS_smaller_correct = (selected_pair1_deltaPhi_muSS_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/deltaPhi_muSS_MET_smaller_correct", deltaPhi_muSS_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test: pick mu_ss with LARGER MT
            bool selected_pair1_MT_muSS_larger = (MT_ss1_MET > MT_ss2_MET);
            bool MT_muSS_larger_correct = (selected_pair1_MT_muSS_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_muSS_MET_larger_correct", MT_muSS_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test: pick mu_ss with SMALLER MT
            bool selected_pair1_MT_muSS_smaller = (MT_ss1_MET < MT_ss2_MET);
            bool MT_muSS_smaller_correct = (selected_pair1_MT_muSS_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_muSS_MET_smaller_correct", MT_muSS_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // MT asymmetry discrimination
            // Calculate MT for mu_os (shared by both pairs)
            float MT_mu_os_MET = sqrt(2 * mu_os.Pt() * METv.Pt() * (1 - cos(mu_os.DeltaPhi(METv))));

            // MT asymmetry for pair1 (mu_ss1 + mu_os)
            float MT_asymmetry1 = abs(MT_ss1_MET - MT_mu_os_MET) / (MT_ss1_MET + MT_mu_os_MET);

            // MT asymmetry for pair2 (mu_ss2 + mu_os)
            float MT_asymmetry2 = abs(MT_ss2_MET - MT_mu_os_MET) / (MT_ss2_MET + MT_mu_os_MET);

            // Test: pick pair with SMALLER MT asymmetry (signal has both muons from A → more similar MT)
            bool selected_pair1_MT_asym_smaller = (MT_asymmetry1 < MT_asymmetry2);
            bool MT_asym_smaller_correct = (selected_pair1_MT_asym_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_asymmetry_smaller_correct", MT_asym_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test: pick pair with LARGER MT asymmetry
            bool selected_pair1_MT_asym_larger = (MT_asymmetry1 > MT_asymmetry2);
            bool MT_asym_larger_correct = (selected_pair1_MT_asym_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/MT_asymmetry_larger_correct", MT_asym_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Invariant mass discrimination
            float mass1 = pair1.M();
            float mass2 = pair2.M();

            // Test: pick pair with SMALLER mass
            bool selected_pair1_mass_smaller = (mass1 < mass2);
            bool mass_smaller_correct = (selected_pair1_mass_smaller == isSignalPair1);
            FillHist(channelStr + "/Discrimination/mass_smaller_correct", mass_smaller_correct ? 1.0 : 0.0, w, 2, 0., 2.);

            // Test: pick pair with LARGER mass
            bool selected_pair1_mass_larger = (mass1 > mass2);
            bool mass_larger_correct = (selected_pair1_mass_larger == isSignalPair1);
            FillHist(channelStr + "/Discrimination/mass_larger_correct", mass_larger_correct ? 1.0 : 0.0, w, 2, 0., 2.);
        }
    }
}
