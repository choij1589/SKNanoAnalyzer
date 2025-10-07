#include "PromptTreeProducer.h"

PromptTreeProducer::PromptTreeProducer() : channel(Channel::NONE) {}

PromptTreeProducer::~PromptTreeProducer() {}

void PromptTreeProducer::WriteHist() {
    GetOutfile()->mkdir("tree");
    GetOutfile()->cd("tree");
    for (auto& [systName, tree] : trees) {
        tree->Write();
    }
    GetOutfile()->cd();
}

void PromptTreeProducer::initializeAnalyzer() {
    TriLeptonBase::initializeAnalyzer();

    // Determine channel
    if (Run1E2Mu) channel = Channel::SR1E2Mu;
    if (Run3Mu) channel = Channel::SR3Mu;
    if (Run1E2Mu && Run3Mu) {
        throw std::runtime_error("Run1E2Mu and Run3Mu cannot be set at the same time");
    }
    if (channel == Channel::NONE) {
        throw std::runtime_error("Run1E2Mu or Run3Mu must be set");
    }

    // Initialize SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", DataStream, DataEra);
    } else if (!RunSyst) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", MCSample, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/TriLeptonSystematics.yaml", MCSample, DataEra);
    }

    // Initialize output trees for each systematic
    GetOutfile()->cd();

    // Collect all systematic names
    std::vector<TString> systNames = {"Central"};

    if (!IsDATA && systHelper) {
        // Add weight-only systematics
        for (const auto& systName : systHelper->getWeightOnlySystematics()) {
            systNames.push_back(systName + "_Up");
            systNames.push_back(systName + "_Down");
        }

        // Add systematics requiring event loop
        for (const auto& syst : *systHelper) {
            if (syst.iter_name != "Central" && systHelper->findSystematic(syst.syst_name)->evtLoopAgain) {
                systNames.push_back(syst.iter_name);
            }
        }
    }

    // Create trees and branches for each systematic
    for (const auto& systName : systNames) {
        TTree* tree = new TTree("Events_" + systName, "");

        tree->Branch("run", &RunNumber);
        tree->Branch("event", &EventNumber);
        tree->Branch("lumi", &LumiBlock);

        tree->Branch("mass1", &mass1[systName]);
        tree->Branch("mass2", &mass2[systName]);

        // GraphNet score branches (placeholder - not functional yet)
        tree->Branch("score_MHc160_MA85_vs_nonprompt", &score_MHc160_MA85_vs_nonprompt[systName]);
        tree->Branch("score_MHc160_MA85_vs_diboson", &score_MHc160_MA85_vs_diboson[systName]);
        tree->Branch("score_MHc160_MA85_vs_ttZ", &score_MHc160_MA85_vs_ttZ[systName]);
        tree->Branch("score_MHc130_MA90_vs_nonprompt", &score_MHc130_MA90_vs_nonprompt[systName]);
        tree->Branch("score_MHc130_MA90_vs_diboson", &score_MHc130_MA90_vs_diboson[systName]);
        tree->Branch("score_MHc130_MA90_vs_ttZ", &score_MHc130_MA90_vs_ttZ[systName]);
        tree->Branch("score_MHc100_MA95_vs_nonprompt", &score_MHc100_MA95_vs_nonprompt[systName]);
        tree->Branch("score_MHc100_MA95_vs_diboson", &score_MHc100_MA95_vs_diboson[systName]);
        tree->Branch("score_MHc100_MA95_vs_ttZ", &score_MHc100_MA95_vs_ttZ[systName]);

        tree->Branch("fold", &fold[systName]);
        tree->Branch("weight", &weight[systName]);

        tree->SetDirectory(0);
        trees[systName] = tree;
    }
}

void PromptTreeProducer::executeEvent() {
    Event ev = GetEvent();

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoVetoMap || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    RVec<Gen> truth = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Initialize tree contents
    initTreeContents();

    // Process events with systematics
    if (!IsDATA && RunSyst && systHelper) {
        // Step 1: Process Central objects and weight-only systematics
        RecoObjects centralObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, "Central");
        Channel selectedChannel = selectEvent(ev, centralObjects, truth, "Central");

        if (selectedChannel != Channel::NONE) {
            // Fill Central with nominal weights
            WeightInfo centralWeights = getWeights(selectedChannel, ev, centralObjects, genJets, "Central");
            fillTree(selectedChannel, centralObjects, centralWeights, "Central", centralObjects.METv);

            // Process weight-only systematics using Central objects
            processWeightOnlySystematics(selectedChannel, ev, centralObjects, genJets, centralObjects.METv);
        }

        // Process systematics requiring evtLoopAgain
        for (const auto& syst : *systHelper) {
            TString systName = syst.iter_name;

            // Skip Central (already processed) and weight-only systematics
            if (systName == "Central" || !systHelper->findSystematic(syst.syst_name)->evtLoopAgain) continue;

            // Define objects with systematic variation
            RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, systName);
            Channel systChannel = selectEvent(ev, recoObjects, truth, systName);

            if (systChannel != Channel::NONE) {
                WeightInfo weights = getWeights(systChannel, ev, recoObjects, genJets, systName);
                // ALWAYS use Central METv for fold calculation
                fillTree(systChannel, recoObjects, weights, systName, centralObjects.METv);
            }
        }
    } else {
        // Process only Central for DATA or when systematics are off
        RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, "Central");
        Channel selectedChannel = selectEvent(ev, recoObjects, truth, "Central");

        if (selectedChannel != Channel::NONE) {
            WeightInfo weights = getWeights(selectedChannel, ev, recoObjects, genJets, "Central");
            fillTree(selectedChannel, recoObjects, weights, "Central", recoObjects.METv);
        }
    }
}

PromptTreeProducer::RecoObjects PromptTreeProducer::defineObjects(Event& ev,
                                                                   const RVec<Muon>& rawMuons,
                                                                   const RVec<Electron>& rawElectrons,
                                                                   const RVec<Jet>& rawJets,
                                                                   const RVec<GenJet>& genJets,
                                                                   const TString& syst) {
    // Create copies of the raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Apply scale variations based on systematic name
    if (syst.Contains("ElectronEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allElectrons = ScaleElectrons(ev, allElectrons, variation);
    } else if (syst.Contains("ElectronRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allElectrons = SmearElectrons(allElectrons, variation);
    } else if (syst.Contains("MuonEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allMuons = ScaleMuons(allMuons, variation);
    } else if (syst.Contains("JetEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allJets = ScaleJets(allJets, variation, "total");
    } else if (syst.Contains("JetRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allJets = SmearJets(allJets, genJets, variation);
    }

    // Get MET and apply Type-I correction
    Particle METv_default;
    if (syst.Contains("UnclusteredEn")) {
        Event::MET_Syst variation = syst.Contains("Up") ? Event::MET_Syst::UE_UP : Event::MET_Syst::UE_DOWN;
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI, variation);
    } else {
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    }
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("tight"), 15., 2.5);

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
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.tightElectrons = tightElectrons;
    objects.jets = jets;
    objects.bjets = bjets;
    objects.METv = METv;

    return objects;
}

PromptTreeProducer::Channel PromptTreeProducer::selectEvent(Event& ev,
                                                             const RecoObjects& recoObjects,
                                                             const RVec<Gen>& truth,
                                                             const TString& syst) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    bool is1E2Mu = (tightElectrons.size() == 1 && vetoElectrons.size() == 1 &&
                    tightMuons.size() == 2 && vetoMuons.size() == 2);
    bool is3Mu = (tightMuons.size() == 3 && vetoMuons.size() == 3 &&
                  tightElectrons.size() == 0 && vetoElectrons.size() == 0);

    if (channel == Channel::SR1E2Mu && !is1E2Mu) return Channel::NONE;
    if (channel == Channel::SR3Mu && !is3Mu) return Channel::NONE;

    // For conversion samples
    if (!IsDATA && (MCSample.Contains("DYJets") || MCSample.Contains("ZGToLLG") || MCSample.Contains("DYGTo2LG"))) {
        RVec<Muon> convMuons;
        RVec<Muon> fakeMuons;
        RVec<Electron> convElectrons;

        for (const auto& mu : tightMuons) {
            int lepType = GetLeptonType(mu, truth);
            if (lepType == 4 || lepType == 5 || lepType == -5 || lepType == -6) convMuons.emplace_back(mu);
            if (lepType == -1 || lepType == -2 || lepType == -3 || lepType == -4) fakeMuons.emplace_back(mu);
        }
        for (const auto& ele : tightElectrons) {
            int lepType = GetLeptonType(ele, truth);
            if (lepType == 4 || lepType == 5 || lepType == -5 || lepType == -6) convElectrons.emplace_back(ele);
        }

        // Remove hadronic contribution
        if (channel == Channel::SR1E2Mu) {
            if (fakeMuons.size() != 0) return Channel::NONE;
            if ((convElectrons.size() + convMuons.size()) == 0) return Channel::NONE;
        }
        if (channel == Channel::SR3Mu) {
            if (fakeMuons.size() != 0) return Channel::NONE;
            if (convMuons.size() == 0) return Channel::NONE;
        }
    }

    // Patching sample
    RVec<Lepton> leptons;
    for (const auto& mu : tightMuons) leptons.emplace_back(mu);
    for (const auto& ele : tightElectrons) leptons.emplace_back(ele);

    bool isLowPT = false;
    for (const auto& l : leptons) {
        if (l.Pt() < 15.) {
            isLowPT = true;
            break;
        }
    }

    if (MCSample.Contains("DYJets") && !isLowPT) return Channel::NONE;
    if ((MCSample.Contains("ZGToLLG") || MCSample.Contains("DYGTo2LG")) && isLowPT) return Channel::NONE;

    // 1E2Mu baseline
    if (channel == Channel::SR1E2Mu) {
        if (!ev.PassTrigger(EMuTriggers)) return Channel::NONE;

        const Muon& mu1 = tightMuons.at(0);
        const Muon& mu2 = tightMuons.at(1);
        const Electron& ele = tightElectrons.at(0);

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

        const Muon& mu1 = tightMuons.at(0);
        const Muon& mu2 = tightMuons.at(1);
        const Muon& mu3 = tightMuons.at(2);

        if (mu1.Pt() <= 20.) return Channel::NONE;
        if (mu2.Pt() <= 10.) return Channel::NONE;
        if (mu3.Pt() <= 10.) return Channel::NONE;
        if (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) != 1) return Channel::NONE;

        auto [pair1, pair2] = makePairs(tightMuons);
        if (pair1.M() <= 12.) return Channel::NONE;
        if (pair2.M() <= 12.) return Channel::NONE;
        if (jets.size() < 2) return Channel::NONE;
        if (bjets.size() < 1) return Channel::NONE;

        return Channel::SR3Mu;
    }

    return Channel::NONE;
}

void PromptTreeProducer::processWeightOnlySystematics(Channel channel,
                                                       const Event& event,
                                                       const RecoObjects& recoObjects,
                                                       const RVec<GenJet>& genJets,
                                                       const Particle& centralMETv) {
    const std::vector<std::string>& allWeightOnlySystematics = systHelper->getWeightOnlySystematics();
    std::vector<std::string> channelFilteredSystematics;

    for (const std::string& systName : allWeightOnlySystematics) {
        bool includeSystematic = true;
        if (systName == "DblMuTrigSF" && channel != Channel::SR3Mu) includeSystematic = false;
        if (systName == "EMuTrigSF" && channel != Channel::SR1E2Mu) includeSystematic = false;
        if (systName == "ElectronIDSF" && channel != Channel::SR1E2Mu) includeSystematic = false;

        if (includeSystematic) {
            channelFilteredSystematics.push_back(systName);
        }
    }

    // Process each weight-only systematic (Up and Down variations)
    // ALWAYS use Central METv for fold calculation
    for (const std::string& systName : channelFilteredSystematics) {
        // Up variation
        TString systNameUp = systName + "_Up";
        WeightInfo weightsUp = getWeights(channel, event, recoObjects, genJets, systNameUp);
        fillTree(channel, recoObjects, weightsUp, systNameUp, centralMETv);

        // Down variation
        TString systNameDown = systName + "_Down";
        WeightInfo weightsDown = getWeights(channel, event, recoObjects, genJets, systNameDown);
        fillTree(channel, recoObjects, weightsDown, systNameDown, centralMETv);
    }
}

PromptTreeProducer::WeightInfo PromptTreeProducer::getWeights(Channel channel,
                                                               const Event& event,
                                                               const RecoObjects& recoObjects,
                                                               const RVec<GenJet>& genJets,
                                                               const TString& syst) {
    WeightInfo weights;

    if (IsDATA) {
        weights.genWeight = 1.;
        weights.prefireWeight = 1.;
        weights.pileupWeight = 1.;
        weights.muonRecoSF = 1.;
        weights.muonIDSF = 1.;
        weights.eleRecoSF = 1.;
        weights.eleIDSF = 1.;
        weights.trigSF = 1.;
        weights.pileupIDSF = 1.;
        weights.btagSF = 1.;
        weights.WZNjetsSF = 1.;
        return weights;
    }

    weights.genWeight = MCweight() * event.GetTriggerLumi("Full");

    MyCorrection::variation var = MyCorrection::variation::nom;
    if (syst.Contains("_Up")) var = MyCorrection::variation::up;
    else if (syst.Contains("_Down")) var = MyCorrection::variation::down;

    // L1 prefire weight
    if (syst.Contains("L1Prefire")) {
        weights.prefireWeight = GetL1PrefireWeight(var);
    } else {
        weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);
    }

    // Pileup weight
    if (syst.Contains("PileupReweight")) {
        weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), var);
    } else {
        weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), MyCorrection::variation::nom);
    }

    // Lepton scale factors
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    const RVec<Muon>& muons = recoObjects.tightMuons;

    weights.muonRecoSF = myCorr->GetMuonRECOSF(muons);

    if (syst.Contains("MuonIDSF")) {
        weights.muonIDSF = myCorr->GetMuonIDSF("TopHNT", muons, var);
    } else {
        weights.muonIDSF = myCorr->GetMuonIDSF("TopHNT", muons);
    }

    weights.eleRecoSF = myCorr->GetElectronRECOSF(electrons);
    if (syst.Contains("ElectronIDSF")) {
        weights.eleIDSF = myCorr->GetElectronIDSF("TopHNT", electrons, var);
    } else {
        weights.eleIDSF = myCorr->GetElectronIDSF("TopHNT", electrons);
    }

    // Trigger scale factors
    if (channel == Channel::SR1E2Mu) {
        if (syst.Contains("EMuTrigSF")) {
            weights.trigSF = myCorr->GetEMuTriggerSF(electrons, muons, var);
        } else {
            weights.trigSF = myCorr->GetEMuTriggerSF(electrons, muons);
        }
    } else if (channel == Channel::SR3Mu) {
        if (syst.Contains("DblMuTrigSF")) {
            weights.trigSF = myCorr->GetDblMuTriggerSF(muons, var);
        } else {
            weights.trigSF = myCorr->GetDblMuTriggerSF(muons);
        }
    } else {
        weights.trigSF = 1.;
    }

    // PileupJet ID scale factor (Run 2 only)
    if (Run == 2) {
        const RVec<Jet>& jets = recoObjects.jets;
        unordered_map<int, int> matched_idx = GenJetMatching(jets, genJets, fixedGridRhoFastjetAll, 0.4, 10.);
        if (syst.Contains("PileupJetIDSF")) {
            weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", var);
        } else {
            weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", MyCorrection::variation::nom);
        }
    } else {
        weights.pileupIDSF = 1.;
    }

    // B-tagging scale factor
    const RVec<Jet>& jets = recoObjects.jets;
    JetTagging::JetFlavTagger tagger = JetTagging::JetFlavTagger::DeepJet;
    JetTagging::JetFlavTaggerWP wp = JetTagging::JetFlavTaggerWP::Medium;
    JetTagging::JetTaggingSFMethod method = JetTagging::JetTaggingSFMethod::mujets;
    string source = "central";
    if (syst.Contains("HFcorr")) source = "hf_corr";
    if (syst.Contains("HFuncorr")) source = "hf_uncorr";
    if (syst.Contains("LFcorr")) source = "lf_corr";
    if (syst.Contains("LFuncorr")) source = "lf_uncorr";
    weights.btagSF = myCorr->GetBTaggingReweightMethod1a(jets, tagger, wp, method, var, source);

    // WZ Njets scale factor
    weights.WZNjetsSF = 1.;
    if (Run == 3 && (MCSample.Contains("WZTo3LNu") || MCSample.Contains("ZZTo4L")) && !RunNoWZSF) {
        float njets = static_cast<float>(jets.size());
        weights.WZNjetsSF = myCorr->GetWZNjetsSF(njets, "Central");
        if (syst.Contains("WZNjetsSF_prompt")) {
            if (syst.Contains("Up")) {
                weights.WZNjetsSF = myCorr->GetWZNjetsSF(njets, "prompt_up");
            } else if (syst.Contains("Down")) {
                weights.WZNjetsSF = myCorr->GetWZNjetsSF(njets, "prompt_down");
            }
        }
        if (syst.Contains("WZNjetsSF_nonprompt")) {
            if (syst.Contains("Up")) {
                weights.WZNjetsSF = myCorr->GetWZNjetsSF(njets, "nonprompt_up");
            } else {
                weights.WZNjetsSF = myCorr->GetWZNjetsSF(njets, "nonprompt_down");
            }
        }
    }

    return weights;
}

void PromptTreeProducer::fillTree(Channel channel, const RecoObjects& recoObjects,
                                   const WeightInfo& weights, const TString& syst,
                                   const Particle& centralMETv) {
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    // Calculate total weight
    if (IsDATA) {
        weight[syst] = 1.;
    } else {
        weight[syst] = weights.genWeight * weights.prefireWeight * weights.pileupWeight *
                       weights.muonRecoSF * weights.muonIDSF * weights.eleRecoSF * weights.eleIDSF *
                       weights.trigSF * weights.pileupIDSF * weights.btagSF * weights.WZNjetsSF;
    }

    // Fill mass variables
    if (channel == Channel::SR1E2Mu) {
        Particle pair = makePair(muons);
        mass1[syst] = pair.M();
        mass2[syst] = -999.;
    } else if (channel == Channel::SR3Mu) {
        auto [pair1, pair2] = makePairs(muons);
        mass1[syst] = pair1.M();
        mass2[syst] = pair2.M();
    }

    // Calculate fold using Central METv (same for all systematics)
    fold[syst] = calculateFold(centralMETv);

    // Evaluate GraphNet scores using the varied METv for this systematic
    const Particle& METv = recoObjects.METv;
    evalScore(muons, electrons, jets, bjets, METv, syst);

    // Fill tree for this systematic
    trees[syst]->Fill();
}

Particle PromptTreeProducer::makePair(const RVec<Muon>& muons) {
    if (muons.size() != 2) {
        throw std::runtime_error("makePair requires exactly 2 muons");
    }
    return muons.at(0) + muons.at(1);
}

std::pair<Particle, Particle> PromptTreeProducer::makePairs(const RVec<Muon>& muons) {
    if (muons.size() != 3) {
        throw std::runtime_error("makePairs requires exactly 3 muons");
    }

    const Muon& mu1 = muons.at(0);
    const Muon& mu2 = muons.at(1);
    const Muon& mu3 = muons.at(2);

    Particle pair1, pair2;
    if (mu1.Charge() == mu2.Charge()) {
        pair1 = mu1 + mu3;
        pair2 = mu2 + mu3;
    } else if (mu1.Charge() == mu3.Charge()) {
        pair1 = mu1 + mu2;
        pair2 = mu3 + mu2;
    } else {  // mu2.Charge() == mu3.Charge()
        pair1 = mu1 + mu2;
        pair2 = mu1 + mu3;
    }

    return {pair1, pair2};
}

void PromptTreeProducer::initTreeContents() {
    for (auto& [systName, tree] : trees) {
        mass1[systName] = -999.;
        mass2[systName] = -999.;
        score_MHc160_MA85_vs_nonprompt[systName] = -999.;
        score_MHc160_MA85_vs_diboson[systName] = -999.;
        score_MHc160_MA85_vs_ttZ[systName] = -999.;
        score_MHc130_MA90_vs_nonprompt[systName] = -999.;
        score_MHc130_MA90_vs_diboson[systName] = -999.;
        score_MHc130_MA90_vs_ttZ[systName] = -999.;
        score_MHc100_MA95_vs_nonprompt[systName] = -999.;
        score_MHc100_MA95_vs_diboson[systName] = -999.;
        score_MHc100_MA95_vs_ttZ[systName] = -999.;
        fold[systName] = -999;
        weight[systName] = -999.;
    }
}

int PromptTreeProducer::calculateFold(const Particle& centralMETv) {
    // Placeholder for fold calculation based on Central METv
    // TODO: Implement when model is trained
    // The fold should be calculated based on int(centralMETv.Pt()) + 1 or similar

    // For now, return dummy value
    return -999;
}

void PromptTreeProducer::evalScore(const RVec<Muon>& muons, const RVec<Electron>& electrons,
                                    const RVec<Jet>& jets, const RVec<Jet>& bjets,
                                    const Particle& METv, const TString& syst) {
    // Placeholder for GraphNet score evaluation using varied METv
    // TODO: Implement when model is trained

    // For now, set dummy values
    score_MHc160_MA85_vs_nonprompt[syst] = -999.;
    score_MHc160_MA85_vs_diboson[syst] = -999.;
    score_MHc160_MA85_vs_ttZ[syst] = -999.;
    score_MHc130_MA90_vs_nonprompt[syst] = -999.;
    score_MHc130_MA90_vs_diboson[syst] = -999.;
    score_MHc130_MA90_vs_ttZ[syst] = -999.;
    score_MHc100_MA95_vs_nonprompt[syst] = -999.;
    score_MHc100_MA95_vs_diboson[syst] = -999.;
    score_MHc100_MA95_vs_ttZ[syst] = -999.;
}
