#include "SimpleDiLepton.h"

SimpleDiLepton::SimpleDiLepton() : channel(Channel::NONE) {}

SimpleDiLepton::~SimpleDiLepton() {}

void SimpleDiLepton::initializeAnalyzer() {
    DiLeptonBase::initializeAnalyzer();

    // Determine channel
    if (RunDiMu) channel = Channel::DIMU;
    if (RunEMu) channel = Channel::EMU;
}

void SimpleDiLepton::executeEvent() {
    Event ev = GetEvent();

    // Initial cutflow entry
    float initialWeight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    fillCutflow(CutStage::Initial, Channel::NONE, initialWeight);

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    fillCutflow(CutStage::NoiseFilter, Channel::NONE, initialWeight);

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoJetVeto || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;
    fillCutflow(CutStage::VetoMap, Channel::NONE, initialWeight);

    RVec<Electron> rawElectrons = GetAllElectrons();

    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Process only Central (no systematics)
    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets);
    fillCutflow(CutStage::LeptonSelection, Channel::NONE, initialWeight);

    Channel selectedChannel = selectEvent(ev, recoObjects);

    if (selectedChannel != Channel::NONE) {
        fillCutflow(CutStage::Final, selectedChannel, initialWeight);
        WeightInfo weights = getWeights(selectedChannel, ev, recoObjects, genParts);
        fillObjects(selectedChannel, recoObjects, weights);
    }
}

SimpleDiLepton::RecoObjects SimpleDiLepton::defineObjects(Event& ev,
                                                         const RVec<Muon>& rawMuons,
                                                         const RVec<Electron>& rawElectrons,
                                                         const RVec<Jet>& rawJets,
                                                         const RVec<GenJet>& genJets) {

    // Create copies of the raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Get MET from event and apply Type-I correction
    Particle METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Define tight leptons with custom ID: POGMediumID + MiniPFRelIso < 0.1
    RVec<Muon> tightMuons;
    for (const auto& mu : allMuons) {
        if (mu.Pt() > 10. && fabs(mu.Eta()) < 2.4 &&
            mu.isPOGMediumId() && mu.MiniPFRelIso() < 0.1) {
            tightMuons.emplace_back(mu);
        }
    }

    RVec<Electron> tightElectrons;
    for (const auto& el : allElectrons) {
        if (el.Pt() > 15. && fabs(el.Eta()) < 2.5 &&
            el.isMVANoIsoWP90() && el.MiniPFRelIso() < 0.1) {
            tightElectrons.emplace_back(el);
        }
    }

    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", 20., max_jeteta);
    if (Run == 2) {
        tightJets = SelectJets(tightJets, "loosePuId", 20., max_jeteta);
        RVec<Jet> tightJets_vetoMap;
        for (const auto &jet: tightJets)
            if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.emplace_back(jet);
        if (!RunNoJetVeto) tightJets = tightJets_vetoMap;
    }
    RVec<Jet> tightJets_vetoLep = JetsVetoLeptonInside(tightJets, tightElectrons, tightMuons, 0.4);

    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (auto& jet : tightJets_vetoLep) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    RecoObjects objects;
    objects.tightMuons = tightMuons;
    objects.tightElectrons = tightElectrons;
    objects.tightJets = tightJets;
    objects.tightJets_vetoLep = tightJets_vetoLep;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv_default = METv_default;
    objects.METv = METv;

    return objects;
}

SimpleDiLepton::Channel SimpleDiLepton::selectEvent(Event& ev, const RecoObjects& recoObjects) {
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;

    float weight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    bool isDiMu = (tightMuons.size() == 2 && tightElectrons.size() == 0);
    bool isEMu = (tightMuons.size() == 1 && tightElectrons.size() == 1);

    if (channel == Channel::DIMU) {
        if (!isDiMu) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::DIMU, weight);
    }
    if (channel == Channel::EMU) {
        if (!isEMu) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::EMU, weight);
    }

    // DiMu selection
    if (channel == Channel::DIMU) {
        if (!ev.PassTrigger(DblMuTriggers)) return Channel::NONE;
        fillCutflow(CutStage::Trigger, Channel::DIMU, weight);

        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        if (! (mu1.Pt() > 20.)) return Channel::NONE;
        if (! (mu2.Pt() > 10.)) return Channel::NONE;
        if (! (mu1.Charge() + mu2.Charge() == 0)) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (! (pair.M() > 50.)) return Channel::NONE;
        fillCutflow(CutStage::KinematicCuts, Channel::DIMU, weight);
        return Channel::DIMU;
    }
    // EMu selection
    else if (channel == Channel::EMU) {
        if (!ev.PassTrigger(EMuTriggers)) return Channel::NONE;
        fillCutflow(CutStage::Trigger, Channel::EMU, weight);

        const Muon& mu = tightMuons[0];
        const Electron& ele = tightElectrons[0];
        if (!((mu.Pt() > 25. && ele.Pt() > 15.) || (mu.Pt() > 10. && ele.Pt() > 25.))) return Channel::NONE;
        if (! (mu.Charge() + ele.Charge() == 0)) return Channel::NONE;
        if (! (mu.DeltaR(ele) > 0.4)) return Channel::NONE;
        fillCutflow(CutStage::KinematicCuts, Channel::EMU, weight);

        if (! (jets.size() >= 2)) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::EMU, weight);

        if (! (bjets.size() >= 1)) return Channel::NONE;
        fillCutflow(CutStage::BjetRequirements, Channel::EMU, weight);

        return Channel::EMU;
    }

    return Channel::NONE;
}

SimpleDiLepton::WeightInfo SimpleDiLepton::getWeights(const SimpleDiLepton::Channel& channel,
                                                     const Event& event,
                                                     const RecoObjects& recoObjects,
                                                     const RVec<Gen>& genParts) {
    WeightInfo weights;

    if (IsDATA) {
        weights.genWeight = 1.;
        weights.prefireWeight = 1.;
        weights.pileupWeight = 1.;
        weights.topPtWeight = 1.;
        return weights;
    }

    // Basic MC weights without ID and trigger scale factors
    weights.genWeight = MCweight() * event.GetTriggerLumi("Full");

    // L1 prefire weight
    weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);

    // Pileup weight
    weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), MyCorrection::variation::nom);

    // Top pT reweighting
    weights.topPtWeight = 1.;
    if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ")) {
        RVec<Gen> genParts = GetAllGens();
        weights.topPtWeight = myCorr->GetTopPtReweight(genParts);
    }

    return weights;
}

void SimpleDiLepton::fillObjects(const SimpleDiLepton::Channel& channel, const RecoObjects& recoObjects,
                                const WeightInfo& weights) {
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle& METv_default = recoObjects.METv_default;
    const Particle& METv = recoObjects.METv;

    TString channelStr = channelToString(channel);
    TString syst = "Central";
    float weight = 1.;

    if (!IsDATA) {
        weight = weights.genWeight * weights.prefireWeight * weights.pileupWeight * weights.topPtWeight;

        FillHist(Form("%s/%s/weights/genWeight", channelStr.Data(), syst.Data()), weights.genWeight, 1., 200, -10000, 10000.);
        FillHist(Form("%s/%s/weights/prefireWeight", channelStr.Data(), syst.Data()), weights.prefireWeight, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/pileupWeight", channelStr.Data(), syst.Data()), weights.pileupWeight, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/topPtWeight", channelStr.Data(), syst.Data()), weights.topPtWeight, 1., 100, -5., 5.);
    }

    // Fill muon histograms
    for (size_t idx = 0; idx < muons.size(); ++idx) {
        const Muon& mu = muons.at(idx);
        FillHist(Form("%s/%s/muons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), mu.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/muons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/muons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/muons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), mu.M(), weight, 10, 0., 1.);
        FillHist(Form("%s/%s/muons/%zu/miniRelIso", channelStr.Data(), syst.Data(), idx+1), mu.MiniPFRelIso(), weight, 100, 0., 1.0);
        FillHist(Form("%s/%s/muons/%zu/sip3d", channelStr.Data(), syst.Data(), idx+1), mu.SIP3D(), weight, 100, 0., 10.);
        FillHist(Form("%s/%s/muons/%zu/dz", channelStr.Data(), syst.Data(), idx+1), mu.dZ(), weight, 100, -0.5, 0.5);
        FillHist(Form("%s/%s/muons/%zu/tkRelIso", channelStr.Data(), syst.Data(), idx+1), mu.TkRelIso(), weight, 100, 0., 1.);
    }

    // Fill electron histograms
    for (size_t idx = 0; idx < electrons.size(); ++idx) {
        const Electron& ele = electrons.at(idx);
        FillHist(Form("%s/%s/electrons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), ele.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/electrons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), ele.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/electrons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/electrons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), ele.M(), weight, 100, 0., 1.);
        FillHist(Form("%s/%s/electrons/%zu/miniRelIso", channelStr.Data(), syst.Data(), idx+1), ele.MiniPFRelIso(), weight, 100, 0., 1.0);
        FillHist(Form("%s/%s/electrons/%zu/sip3d", channelStr.Data(), syst.Data(), idx+1), ele.SIP3D(), weight, 100, 0., 10.);
        FillHist(Form("%s/%s/electrons/%zu/dz", channelStr.Data(), syst.Data(), idx+1), ele.dZ(), weight, 100, -0.5, 0.5);
        FillHist(Form("%s/%s/electrons/%zu/tkRelIso", channelStr.Data(), syst.Data(), idx+1), ele.TkRelIso(), weight, 100, 0., 1.);
    }

    // Fill jet histograms
    for (size_t idx = 0; idx < jets.size(); ++idx) {
        const Jet& jet = jets.at(idx);
        FillHist(Form("%s/%s/jets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), jet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/rawPt", channelStr.Data(), syst.Data(), idx+1), jet.GetRawPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/originalPt", channelStr.Data(), syst.Data(), idx+1), jet.GetOriginalPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), jet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/jets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), jet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/jets/%zu/mass", channelStr.Data(), syst.Data(), idx+1), jet.M(), weight, 100, 0., 100.);
    }

    // Fill bjet histograms
    for (size_t idx = 0; idx < bjets.size(); ++idx) {
        const Jet& bjet = bjets.at(idx);
        FillHist(Form("%s/%s/bjets/%zu/pt", channelStr.Data(), syst.Data(), idx+1), bjet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/bjets/%zu/eta", channelStr.Data(), syst.Data(), idx+1), bjet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/bjets/%zu/phi", channelStr.Data(), syst.Data(), idx+1), bjet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/bjets/%zu/mass", channelStr.Data(), syst.Data(), idx+1), bjet.M(), weight, 100, 0., 100.);
    }

    FillHist(Form("%s/%s/jets/size", channelStr.Data(), syst.Data()), jets.size(), weight, 20, 0., 20.);
    FillHist(Form("%s/%s/bjets/size", channelStr.Data(), syst.Data()), bjets.size(), weight, 15, 0., 15.);
    FillHist(Form("%s/%s/METv/pt", channelStr.Data(), syst.Data()), METv.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv/phi", channelStr.Data(), syst.Data()), METv.Phi(), weight, 64, -3.2, 3.2);
    FillHist(Form("%s/%s/METv_default/pt", channelStr.Data(), syst.Data()), METv_default.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv_default/phi", channelStr.Data(), syst.Data()), METv_default.Phi(), weight, 64, -3.2, 3.2);

    if (channel == Channel::DIMU) {
        Particle pair = muons.at(0) + muons.at(1);
        FillHist(Form("%s/%s/pair/pt", channelStr.Data(), syst.Data()), pair.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/pair/eta", channelStr.Data(), syst.Data()), pair.Eta(), weight, 100, -5., 5.);
        FillHist(Form("%s/%s/pair/phi", channelStr.Data(), syst.Data()), pair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/pair/mass", channelStr.Data(), syst.Data()), pair.M(), weight, 300, 0., 300.);
    }

    if (! (jets.size() > 2)) return;
    TString subChannelStr = channelStr+"_2j";
    // Fill muon histograms
    for (size_t idx = 0; idx < muons.size(); ++idx) {
        const Muon& mu = muons.at(idx);
        FillHist(Form("%s/%s/muons/%zu/pt", subChannelStr.Data(), syst.Data(), idx+1), mu.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/muons/%zu/eta", subChannelStr.Data(), syst.Data(), idx+1), mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/muons/%zu/phi", subChannelStr.Data(), syst.Data(), idx+1), mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/muons/%zu/mass", subChannelStr.Data(), syst.Data(), idx+1), mu.M(), weight, 10, 0., 1.);
    }

    // Fill electron histograms
    for (size_t idx = 0; idx < electrons.size(); ++idx) {
        const Electron& ele = electrons.at(idx);
        FillHist(Form("%s/%s/electrons/%zu/pt", subChannelStr.Data(), syst.Data(), idx+1), ele.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/electrons/%zu/eta", subChannelStr.Data(), syst.Data(), idx+1), ele.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/electrons/%zu/phi", subChannelStr.Data(), syst.Data(), idx+1), ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/electrons/%zu/mass", subChannelStr.Data(), syst.Data(), idx+1), ele.M(), weight, 100, 0., 1.);
    }

    // Fill jet histograms
    for (size_t idx = 0; idx < jets.size(); ++idx) {
        const Jet& jet = jets.at(idx);
        FillHist(Form("%s/%s/jets/%zu/pt", subChannelStr.Data(), syst.Data(), idx+1), jet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/rawPt", subChannelStr.Data(), syst.Data(), idx+1), jet.GetRawPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/originalPt", subChannelStr.Data(), syst.Data(), idx+1), jet.GetOriginalPt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/jets/%zu/eta", subChannelStr.Data(), syst.Data(), idx+1), jet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/jets/%zu/phi", subChannelStr.Data(), syst.Data(), idx+1), jet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/jets/%zu/mass", subChannelStr.Data(), syst.Data(), idx+1), jet.M(), weight, 100, 0., 100.);
    }

    // Fill bjet histograms
    for (size_t idx = 0; idx < bjets.size(); ++idx) {
        const Jet& bjet = bjets.at(idx);
        FillHist(Form("%s/%s/bjets/%zu/pt", subChannelStr.Data(), syst.Data(), idx+1), bjet.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/bjets/%zu/eta", subChannelStr.Data(), syst.Data(), idx+1), bjet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/bjets/%zu/phi", subChannelStr.Data(), syst.Data(), idx+1), bjet.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/bjets/%zu/mass", subChannelStr.Data(), syst.Data(), idx+1), bjet.M(), weight, 100, 0., 100.);
    }

    FillHist(Form("%s/%s/jets/size", subChannelStr.Data(), syst.Data()), jets.size(), weight, 20, 0., 20.);
    FillHist(Form("%s/%s/bjets/size", subChannelStr.Data(), syst.Data()), bjets.size(), weight, 15, 0., 15.);
    FillHist(Form("%s/%s/METv/pt", subChannelStr.Data(), syst.Data()), METv.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv/phi", subChannelStr.Data(), syst.Data()), METv.Phi(), weight, 64, -3.2, 3.2);
    FillHist(Form("%s/%s/METv_default/pt", subChannelStr.Data(), syst.Data()), METv_default.Pt(), weight, 300, 0., 300.);
    FillHist(Form("%s/%s/METv_default/phi", subChannelStr.Data(), syst.Data()), METv_default.Phi(), weight, 64, -3.2, 3.2);

    if (channel == Channel::DIMU) {
        Particle pair = muons.at(0) + muons.at(1);
        FillHist(Form("%s/%s/pair/pt", subChannelStr.Data(), syst.Data()), pair.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/pair/eta", subChannelStr.Data(), syst.Data()), pair.Eta(), weight, 100, -5., 5.);
        FillHist(Form("%s/%s/pair/phi", subChannelStr.Data(), syst.Data()), pair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/pair/mass", subChannelStr.Data(), syst.Data()), pair.M(), weight, 300, 0., 300.);
    }
}

void SimpleDiLepton::fillCutflow(CutStage stage, const Channel& channel, float weight) {
    TString channelStr = channelToString(channel);
    if (channelStr == "NONE") channelStr = "ALL";
    TString syst = "Central";

    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/cutflow", channelStr.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}
