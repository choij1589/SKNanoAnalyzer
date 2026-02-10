#include "MeasFakeRateV4.h"

MeasFakeRateV4::MeasFakeRateV4() : leptonType(LeptonType::NONE), lowestPtCut(10.) {}

MeasFakeRateV4::~MeasFakeRateV4() {}

void MeasFakeRateV4::initializeAnalyzer() {
    // Set flags
    MeasFakeMu = HasFlag("MeasFakeMu");
    MeasFakeEl = HasFlag("MeasFakeEl");
    RunSyst = HasFlag("RunSyst");
    RunNoHEMVeto = HasFlag("RunNoHEMVeto");

    // Determine lepton type and configure triggers
    if (MeasFakeMu) {
        leptonType = LeptonType::MUON;
        ptcorr_bins = {10., 12., 14., 17., 20., 30., 50., 100., 200.};
        abseta_bins = {0., 0.9, 1.6, 2.4};
        triggers = {
            {"HLT_Mu8_TrkIsoVVL",  "Mu8",  10., false},
            {"HLT_Mu17_TrkIsoVVL", "Mu17", 20., false}
        };
        lowestPtCut = 10.;
    } else if (MeasFakeEl) {
        leptonType = LeptonType::ELECTRON;
        ptcorr_bins = {15., 17., 20., 25., 35., 50., 100., 200.};
        abseta_bins = {0., 0.8, 1.479, 2.5};
        triggers = {
            {"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",  "El8",  10., false},
            {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", "El12", 15., false},
            {"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", "El23", 25., false}
        };
        lowestPtCut = 10.;
    } else {
        throw std::runtime_error("[MeasFakeRateV4::initializeAnalyzer] No lepton type specified by flags. Use MeasFakeMu or MeasFakeEl.");
    }

    // Set IDs
    MuonIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    ElectronIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));

    // Initialize correction
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    // Initialize SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/FakeSystematics.data.yaml", DataStream, DataEra);
    } else {
        if (RunSyst) {
            systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/FakeSystematics.mc.yaml", MCSample, DataEra);
        } else {
            systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", DataStream, DataEra);
        }
    }
}

void MeasFakeRateV4::executeEvent() {
    Event ev = GetEvent();

    // Initial cutflow entry
    float initialWeight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    fillCutflow(CutStage::Initial, Channel::NONE, "event", initialWeight, "Central");

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    fillCutflow(CutStage::NoiseFilter, Channel::NONE, "event", initialWeight, "Central");

    RVec<Muon> rawMuons = GetAllMuons();
    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;
    fillCutflow(CutStage::VetoMap, Channel::NONE, "event", initialWeight, "Central");

    RVec<Electron> rawElectrons = GetAllElectrons();

    // Check which triggers fired (store in triggers[i].fired)
    bool anyFired = false;
    for (auto& trig : triggers) {
        trig.fired = ev.PassTrigger(trig.name);
        anyFired |= trig.fired;
    }
    if (!anyFired) return;
    fillCutflow(CutStage::AnyTrigger, Channel::NONE, "event", initialWeight, "Central");

    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Loop over IDs (loose and tight)
    RVec<TString> IDs = {"loose", "tight"};

    for (const auto& ID : IDs) {
        if (RunSyst && systHelper) {
            // Process Central objects and weight-only systematics
            RecoObjects centralObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, centralObjects, ID, "Central");

            if (selectedChannel != Channel::NONE) {
                fillCutflow(CutStage::Final, selectedChannel, ID, initialWeight, "Central");
                WeightInfo centralWeights = getWeights(selectedChannel, ID, ev, centralObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, centralObjects, centralWeights, "Central");

                // process weight-only systematics with central objects
                vector<string> weightOnlySysts = systHelper->getWeightOnlySystematics();
                for (const auto &systName : weightOnlySysts) {
                    TString systNameUp = systName + "_Up";
                    WeightInfo weightsUp = getWeights(selectedChannel, ID, ev, centralObjects, genParts, systNameUp);
                    fillObjects(selectedChannel, ID, centralObjects, weightsUp, systNameUp);

                    TString systNameDown = systName + "_Down";
                    WeightInfo weightsDown = getWeights(selectedChannel, ID, ev, centralObjects, genParts, systNameDown);
                    fillObjects(selectedChannel, ID, centralObjects, weightsDown, systNameDown);
                }
            }

            // Process systematics requiring evtLoopAgain
            for (const auto& syst : *systHelper) {
                TString systName = syst.iter_name;

                // Skip Central and weight-only systematics
                if (systName == "Central" || (!systHelper->findSystematic(syst.syst_name)->evtLoopAgain)) continue;

                RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, systName);
                Channel systChannel = selectEvent(ev, recoObjects, ID, systName);

                if (systChannel != Channel::NONE) {
                    WeightInfo weights = getWeights(systChannel, ID, ev, recoObjects, genParts, systName);
                    fillObjects(systChannel, ID, recoObjects, weights, systName);
                }
            }
        } else {
            // systematics are off
            RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, recoObjects, ID, "Central");

            if (selectedChannel != Channel::NONE) {
                fillCutflow(CutStage::Final, selectedChannel, ID, initialWeight, "Central");
                WeightInfo weights = getWeights(selectedChannel, ID, ev, recoObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, recoObjects, weights, "Central");
            }
        }
    }
}

MeasFakeRateV4::RecoObjects MeasFakeRateV4::defineObjects(Event& ev,
                                                     const RVec<Muon>& rawMuons,
                                                     const RVec<Electron>& rawElectrons,
                                                     const RVec<Jet>& rawJets,
                                                     const RVec<GenJet>& genJets,
                                                     const TString& ID,
                                                     const TString& syst) {
    // Create copies for systematic variations
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Apply systematic variations
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
        TString systSource = "total";
        allJets = ScaleJets(allJets, variation, systSource);
    } else if (syst.Contains("JetRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allJets = SmearJets(allJets, genJets, variation);
    } else {
        // No scale variation
    }

    // Get MET with unclustered energy variation if applicable
    Particle METv_default;
    if (syst.Contains("UnclusteredEn")) {
        Event::MET_Syst variation = syst.Contains("Up") ? Event::MET_Syst::UE_UP : Event::MET_Syst::UE_DOWN;
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI, variation);
    } else {
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    }
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort objects by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select objects based on ID
    // HEM veto: apply to measurement electrons in 2018 unless RunNoHEMVeto flag is set
    bool applyHEMVeto = DataEra.Contains("2018") && !RunNoHEMVeto;

    // Veto leptons: 10 GeV, loose ID (for event selection veto, no HEM veto)
    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);

    // Measurement leptons: muons at 10 GeV, electrons at 15 GeV (with HEM veto for 2018)
    RVec<Muon> looseMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Electron> looseElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("loose"), 15., 2.5, applyHEMVeto);

    RVec<Muon> tightMuons;
    RVec<Electron> tightElectrons;

    if (ID == "loose") {
        tightMuons = looseMuons;
        tightElectrons = looseElectrons;
    } else if (ID == "tight") {
        tightMuons = SelectMuons(allMuons, MuonIDs->GetID("tight"), 10., 2.4);
        tightElectrons = SelectElectrons(looseElectrons, ElectronIDs->GetID("tight"), 15., 2.5, applyHEMVeto);
    }

    // Jet selection with potential selection variations
    // Order: tight -> vetoLep -> vetoMap -> PUID
    const float jetPtCut = getJetPtCut(syst);
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", jetPtCut, max_jeteta);
    tightJets = JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4);
    RVec<Jet> tightJets_noPUID;  // Jets before PUID for SF calculation
    if (Run == 2) {
        RVec<Jet> tightJets_vetoMap;
        for (const auto &jet : tightJets) {
            if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.push_back(jet);
        }
        // Store jets before PUID for SF calculation, then apply PUID
        tightJets_noPUID = tightJets_vetoMap;
        tightJets = SelectJets(tightJets_vetoMap, "loosePuId", jetPtCut, max_jeteta);
    }

    // B-jet selection
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& jet : tightJets) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    // Compute mother jet flavours for leptons (MC only)
    RVec<int> looseMuonJetFlavours, tightMuonJetFlavours;
    RVec<int> looseElectronJetFlavours, tightElectronJetFlavours;

    if (!IsDATA) {
        for (const auto& mu : looseMuons) {
            looseMuonJetFlavours.push_back(getMotherJetFlavour(mu, allJets));
        }
        for (const auto& mu : tightMuons) {
            tightMuonJetFlavours.push_back(getMotherJetFlavour(mu, allJets));
        }
        for (const auto& el : looseElectrons) {
            looseElectronJetFlavours.push_back(getMotherJetFlavour(el, allJets));
        }
        for (const auto& el : tightElectrons) {
            tightElectronJetFlavours.push_back(getMotherJetFlavour(el, allJets));
        }
    }

    RecoObjects objects;
    objects.vetoMuons = vetoMuons;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.looseMuonJetFlavours = looseMuonJetFlavours;
    objects.tightMuonJetFlavours = tightMuonJetFlavours;
    objects.looseElectronJetFlavours = looseElectronJetFlavours;
    objects.tightElectronJetFlavours = tightElectronJetFlavours;
    objects.tightJets = tightJets;
    objects.tightJets_noPUID = tightJets_noPUID;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv = METv;

    return objects;
}

MeasFakeRateV4::Channel MeasFakeRateV4::selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst) {
    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    float weight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");

    if (leptonType == LeptonType::MUON) {
        const bool sglMu = (muons.size() == 1 && vetoMuons.size() == 1 &&
                            electrons.size() == 0 && vetoElectrons.size() == 0);
        const bool dblMu = (muons.size() == 2 && vetoMuons.size() == 2 &&
                            electrons.size() == 0 && vetoElectrons.size() == 0);

        if (! (sglMu || dblMu)) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::INCLUSIVE, ID, weight, syst);

        // Use loosest pT cut for initial selection (trigger-specific cut applied in fillObjects)
        if (! (muons[0].Pt() > lowestPtCut)) return Channel::NONE;
        if (! (jets.size() > 0)) return Channel::NONE;
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::INCLUSIVE, ID, weight, syst);

        if (sglMu) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(muons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            fillCutflow(CutStage::AwayJetRequirements, Channel::INCLUSIVE, ID, weight, syst);
            return Channel::INCLUSIVE;
        } else { // doubleMu
            Particle ZCand = muons[0] + muons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            fillCutflow(CutStage::ZMassWindow, Channel::ZENRICHED, ID, weight, syst);
            return Channel::ZENRICHED;
        }
    } else { // ELECTRON
        const bool sglEl = (electrons.size() == 1 && vetoElectrons.size() == 1 &&
                        muons.size() == 0 && vetoMuons.size() == 0);
        const bool dblEl = (electrons.size() == 2 && vetoElectrons.size() == 2 &&
                        muons.size() == 0 && vetoMuons.size() == 0);

        if (! (sglEl || dblEl)) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::INCLUSIVE, ID, weight, syst);

        // Use loosest pT cut for initial selection (trigger-specific cut applied in fillObjects)
        if (! (electrons[0].Pt() > lowestPtCut)) return Channel::NONE;
        if (! (jets.size() > 0)) return Channel::NONE;
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::INCLUSIVE, ID, weight, syst);

        if (sglEl) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(electrons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            fillCutflow(CutStage::AwayJetRequirements, Channel::INCLUSIVE, ID, weight, syst);
            return Channel::INCLUSIVE;
        } else { // dblEl
            Particle ZCand = electrons[0] + electrons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            fillCutflow(CutStage::ZMassWindow, Channel::ZENRICHED, ID, weight, syst);
            return Channel::ZENRICHED;
        }
    }

    return Channel::NONE;
}

MeasFakeRateV4::WeightInfo MeasFakeRateV4::getWeights(const Channel& channel,
                                                  const TString& ID,
                                                  const Event& event,
                                                  const RecoObjects& recoObjects,
                                                  const RVec<Gen>& genParts,
                                                  const TString& syst) {
    WeightInfo weights;
    weights.genWeight = 1.0;
    weights.prefireWeight = 1.0;
    weights.pileupWeight = 1.0;
    weights.topPtWeight = 1.0;
    weights.muonRecoSF = 1.0;
    weights.eleRecoSF = 1.0;
    weights.btagSF = 1.0;
    weights.pileupIDSF = 1.0;

    if (!IsDATA) {
        weights.genWeight = MCweight() * event.GetTriggerLumi("Full");

        // Determine systematic variation
        MyCorrection::variation var = MyCorrection::variation::nom;
        if (syst.Contains("_Up")) var = MyCorrection::variation::up;
        else if (syst.Contains("_Down")) var = MyCorrection::variation::down;

        // L1 Prefire
        if (syst.Contains("L1Prefire")) {
            weights.prefireWeight = GetL1PrefireWeight(var);
        } else {
            weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);
        }

        // Top pT reweight
        if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ")) {
            weights.topPtWeight = myCorr->GetTopPtReweight(genParts);
        }

        // Lepton reconstruction SF
        const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
        const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;

        if (syst.Contains("MuonRecoSF")) {
            weights.muonRecoSF = myCorr->GetMuonRECOSF(muons, var);
        } else {
            weights.muonRecoSF = myCorr->GetMuonRECOSF(muons, MyCorrection::variation::nom);
        }

        if (syst.Contains("ElectronRecoSF")) {
            weights.eleRecoSF = myCorr->GetElectronRECOSF(electrons, var);
        } else {
            weights.eleRecoSF = myCorr->GetElectronRECOSF(electrons, MyCorrection::variation::nom);
        }

        // B-tagging SF for RequireHeavyTag
        if (syst.Contains("RequireHeavyTag")) {
            weights.btagSF = myCorr->GetBTaggingReweightMethod1a(recoObjects.tightJets,
                                                               JetTagging::JetFlavTagger::DeepJet,
                                                               JetTagging::JetFlavTaggerWP::Medium,
                                                               JetTagging::JetTaggingSFMethod::mujets,
                                                               MyCorrection::variation::nom);
        }

        // Jet PUID SF for Run2
        // Use jets before PUID selection for SF calculation
        if (Run == 2) {
            const RVec<Jet>& jets = recoObjects.tightJets_noPUID;
            const RVec<GenJet>& genJets = recoObjects.genJets;

            unordered_map<int, int> matched_idx = GenJetMatching(jets, genJets, fixedGridRhoFastjetAll, 0.4, 10.);
            if (syst.Contains("PileupJetIDSF")) {
                weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", var);
            } else {
                weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", MyCorrection::variation::nom);
            }
        }

        // Pileup reweighting
        if (syst.Contains("PileupReweight")) {
            weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), var);
        } else {
            weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), MyCorrection::variation::nom);
        }
    }

    return weights;
}

void MeasFakeRateV4::fillObjects(const Channel& channel,
                               const TString& ID,
                               const RecoObjects& recoObjects,
                               const WeightInfo& weights,
                               const TString& syst) {
    float totalWeight = 1.;
    if (!IsDATA) {
        totalWeight = weights.genWeight;
        totalWeight *= weights.prefireWeight;
        totalWeight *= weights.pileupWeight;
        totalWeight *= weights.topPtWeight;
        totalWeight *= weights.muonRecoSF;
        totalWeight *= weights.eleRecoSF;
        if (syst == "RequireHevayTag") totalWeight *= weights.btagSF;
        if (Run == 2) totalWeight *= weights.pileupIDSF;
    }

    // Memory optimization: For promptNormOnly systematics, only fill ZEnriched Z mass
    // This reduces histogram count for weight-only systematic variations
    // Full measurement systematics (MotherJetPt, RequireHeavyTag, etc.) need all histograms
    if (syst != "Central") {
        // Check if this systematic is promptNormOnly
        bool isPromptNormOnly = false;
        if (systHelper) {
            // Extract base syst name (remove _Up/_Down suffix)
            TString baseSyst = syst;
            baseSyst.ReplaceAll("_Up", "").ReplaceAll("_Down", "");
            SystematicHelper::SYST* systInfo = systHelper->findSystematic(baseSyst.Data());
            if (systInfo) {
                isPromptNormOnly = systInfo->promptNormOnly;
            }
        }

        if (isPromptNormOnly) {
            if (channel != Channel::ZENRICHED) return;  // Skip Inclusive channel

            // Only fill Z mass histogram for normalization
            for (const auto& trig : triggers) {
                if (!trig.fired) continue;

                TString basePrefix = channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;

                if (leptonType == LeptonType::MUON) {
                    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
                    float ptcorr_lead = muons[0].Pt()*(1.+max(0., muons[0].MiniPFRelIso()-0.1));
                    if (muons[0].Pt() < trig.trigSafePtCut) continue;

                    const Particle ZCand = muons[0] + muons[1];
                    FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                } else { // ELECTRON
                    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
                    float ptcorr_lead = electrons[0].Pt()*(1.+max(0., electrons[0].MiniPFRelIso()-0.1));
                    if (electrons[0].Pt() < trig.trigSafePtCut) continue;

                    const Particle ZCand = electrons[0] + electrons[1];
                    FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                }
            }
            return;  // Skip all other histograms for promptNormOnly systematics
        }
        // For non-promptNormOnly systematics, fall through to full histogram filling
    }

    // Central systematic: fill all histograms
    const RVec<Jet> &jets = recoObjects.tightJets;
    const RVec<Jet> &bjets = recoObjects.bjets;
    const Particle METv = recoObjects.METv;

    // Loop over triggers that fired and fill histograms for each
    for (const auto& trig : triggers) {
        if (!trig.fired) continue;

        // Build base prefix: Channel/TrigPrefix/ID/Syst
        TString basePrefix = channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;

        // Fill histograms based on channel and lepton type
        if (channel == Channel::INCLUSIVE) {
            if (leptonType == LeptonType::MUON) {
                const Muon& mu = (ID == "loose") ? recoObjects.looseMuons[0] : recoObjects.tightMuons[0];
                float ptcorr = mu.Pt()*(1.+max(0., mu.MiniPFRelIso()-0.1));

                // Check if lepton passes this trigger's pT threshold
                if (mu.Pt() < trig.trigSafePtCut) continue;

                float abseta = fabs(mu.Eta());
                float mT = TMath::Sqrt(2.*mu.Pt()*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
                float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
                TString binName = getBinPrefix(ptcorr, abseta);

                // Fill inclusive channel histograms
                FillHist(basePrefix + "/muon/pt", mu.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muon/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muon/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(basePrefix + "/muon/abseta", abseta, totalWeight, abseta_bins);
                FillHist(basePrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms in Inclusive/ (MC only)
                if (!IsDATA) {
                    int jetFlavour = (ID == "loose") ? recoObjects.looseMuonJetFlavours[0] : recoObjects.tightMuonJetFlavours[0];
                    TString flavourSubdir = getFlavourSubdir(jetFlavour);
                    TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                    FillHist(flavourPrefix + "/muon/pt", mu.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(flavourPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                    FillHist(flavourPrefix + "/muon/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
                    FillHist(flavourPrefix + "/muon/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                    FillHist(flavourPrefix + "/muon/abseta", abseta, totalWeight, abseta_bins);
                    FillHist(flavourPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                    FillHist(flavourPrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                    FillHist(flavourPrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                }

                // Fill binned histograms: binName/Channel/TrigPrefix/ID/Syst
                TString binnedPrefix = binName + "/" + channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;
                FillHist(binnedPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                FillHist(binnedPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                FillHist(binnedPrefix + "/MT", mT, totalWeight, 500, 0., 500.);
                FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 600, 0., 300.);

                // Fill subchannel - determine QCDEnriched or WEnriched
                TString subchannel = "";
                if (mT < 25. && METv.Pt() < 25.) {
                    subchannel = "QCDEnriched";
                } else if (mT > 60.) {
                    subchannel = "WEnriched";
                }

                if (subchannel != "") {
                    TString subchannelPrefix = binName + "/" + subchannel + "/" + trig.prefix + "/" + ID + "/" + syst;
                    FillHist(subchannelPrefix + "/muon/pt", mu.Pt(), totalWeight, 200, 0., 200.);
                    FillHist(subchannelPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                    FillHist(subchannelPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                    FillHist(subchannelPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                    FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);

                    // Fill flavour-binned histograms (MC only)
                    if (!IsDATA) {
                        int jetFlavour = (ID == "loose") ? recoObjects.looseMuonJetFlavours[0] : recoObjects.tightMuonJetFlavours[0];
                        TString flavourSubdir = getFlavourSubdir(jetFlavour);

                        // Flavour-binned histograms under binnedPrefix
                        TString flavourBinnedPrefix = binnedPrefix + "/" + flavourSubdir;
                        FillHist(flavourBinnedPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                        FillHist(flavourBinnedPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                        FillHist(flavourBinnedPrefix + "/MT", mT, totalWeight, 500, 0., 500.);
                        FillHist(flavourBinnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MET", METv.Pt(), totalWeight, 600, 0., 300.);

                        // Flavour-binned histograms under subchannelPrefix
                        TString flavourSubchPrefix = subchannelPrefix + "/" + flavourSubdir;
                        FillHist(flavourSubchPrefix + "/muon/pt", mu.Pt(), totalWeight, 200, 0., 200.);
                        FillHist(flavourSubchPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                        FillHist(flavourSubchPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                        FillHist(flavourSubchPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                        FillHist(flavourSubchPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
                    }
                }
            } else { // ELECTRON
                const Electron& el = (ID=="loose") ? recoObjects.looseElectrons[0]: recoObjects.tightElectrons[0];
                float ptcorr = el.Pt()*(1.+max(0., el.MiniPFRelIso()-0.1));

                // Check if lepton passes this trigger's pT threshold
                if (el.Pt() < trig.trigSafePtCut) continue;

                float abseta = fabs(el.scEta());
                float mT = TMath::Sqrt(2.*el.Pt()*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
                float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
                TString binName = getBinPrefix(ptcorr, abseta);

                // Fill inclusive channel histograms
                FillHist(basePrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electron/phi", el.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(basePrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                FillHist(basePrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms in Inclusive/ (MC only)
                if (!IsDATA) {
                    int jetFlavour = (ID == "loose") ? recoObjects.looseElectronJetFlavours[0] : recoObjects.tightElectronJetFlavours[0];
                    TString flavourSubdir = getFlavourSubdir(jetFlavour);
                    TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                    FillHist(flavourPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(flavourPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                    FillHist(flavourPrefix + "/electron/phi", el.Phi(), totalWeight, 64, -3.2, 3.2);
                    FillHist(flavourPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                    FillHist(flavourPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                    FillHist(flavourPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                    FillHist(flavourPrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                    FillHist(flavourPrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                }

                // Fill binned histograms: binName/Channel/TrigPrefix/ID/Syst
                TString binnedPrefix = binName + "/" + channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;
                FillHist(binnedPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(binnedPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                FillHist(binnedPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);

                // Fill subchannel - determine QCDEnriched or WEnriched
                TString subchannel = "";
                if (mT < 25. && METv.Pt() < 25.) {
                    subchannel = "QCDEnriched";
                } else if (mT > 60.) {
                    subchannel = "WEnriched";
                }

                if (subchannel != "") {
                    TString subchannelPrefix = binName + "/" + subchannel + "/" + trig.prefix + "/" + ID + "/" + syst;
                    FillHist(subchannelPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                    FillHist(subchannelPrefix + "/electron/ptcorr", ptcorr, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/electron/abseta", abseta, totalWeight, 25, 0., 2.5);
                    FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);

                    // Fill flavour-binned histograms (MC only)
                    if (!IsDATA) {
                        int jetFlavour = (ID == "loose") ? recoObjects.looseElectronJetFlavours[0] : recoObjects.tightElectronJetFlavours[0];
                        TString flavourSubdir = getFlavourSubdir(jetFlavour);

                        // Flavour-binned histograms under binnedPrefix
                        TString flavourBinnedPrefix = binnedPrefix + "/" + flavourSubdir;
                        FillHist(flavourBinnedPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                        FillHist(flavourBinnedPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                        FillHist(flavourBinnedPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);

                        // Flavour-binned histograms under subchannelPrefix
                        TString flavourSubchPrefix = subchannelPrefix + "/" + flavourSubdir;
                        FillHist(flavourSubchPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                        FillHist(flavourSubchPrefix + "/electron/ptcorr", ptcorr, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/electron/abseta", abseta, totalWeight, 25, 0., 2.5);
                        FillHist(flavourSubchPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
                    }
                }
            }
        } else if (channel == Channel::ZENRICHED) {
            if (leptonType == LeptonType::MUON) {
                const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;

                // For ZENRICHED, check if leading muon passes trigger pT threshold
                float ptcorr_lead = muons[0].Pt()*(1.+max(0., muons[0].MiniPFRelIso()-0.1));
                if (muons[0].Pt() < trig.trigSafePtCut) continue;

                const Particle ZCand = muons[0] + muons[1];
                FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                FillHist(basePrefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
                FillHist(basePrefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muons/1/pt", muons[0].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muons/1/eta", muons[0].Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muons/1/phi", muons[0].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muons/2/pt", muons[1].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muons/2/eta", muons[1].Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muons/2/phi", muons[1].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms for ZENRICHED (MC only)
                if (!IsDATA) {
                    const RVec<int>& jetFlavours = (ID == "loose") ? recoObjects.looseMuonJetFlavours : recoObjects.tightMuonJetFlavours;
                    for (size_t i = 0; i < muons.size() && i < jetFlavours.size(); i++) {
                        TString flavourSubdir = getFlavourSubdir(jetFlavours[i]);
                        TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                        FillHist(flavourPrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                        FillHist(flavourPrefix + Form("/muons/%zu/pt", i+1), muons[i].Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourPrefix + Form("/muons/%zu/eta", i+1), muons[i].Eta(), totalWeight, 48, -2.4, 2.4);
                    }
                }
            } else { // ELECTRON
                const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;

                // For ZENRICHED, check if leading electron passes trigger pT threshold
                float ptcorr_lead = electrons[0].Pt()*(1.+max(0., electrons[0].MiniPFRelIso()-0.1));
                if (electrons[0].Pt() < trig.trigSafePtCut) continue;

                const Particle ZCand = electrons[0] + electrons[1];
                FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                FillHist(basePrefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
                FillHist(basePrefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electrons/1/pt", electrons[0].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electrons/1/scEta", electrons[0].scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electrons/1/phi", electrons[0].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electrons/2/pt", electrons[1].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electrons/2/scEta", electrons[1].scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electrons/2/phi", electrons[1].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms for ZENRICHED (MC only)
                if (!IsDATA) {
                    const RVec<int>& jetFlavours = (ID == "loose") ? recoObjects.looseElectronJetFlavours : recoObjects.tightElectronJetFlavours;
                    for (size_t i = 0; i < electrons.size() && i < jetFlavours.size(); i++) {
                        TString flavourSubdir = getFlavourSubdir(jetFlavours[i]);
                        TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                        FillHist(flavourPrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                        FillHist(flavourPrefix + Form("/electrons/%zu/pt", i+1), electrons[i].Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourPrefix + Form("/electrons/%zu/scEta", i+1), electrons[i].scEta(), totalWeight, 50, -2.5, 2.5);
                    }
                }
            }
        }
    }  // end trigger loop
}

TString MeasFakeRateV4::getBinPrefix(const double ptcorr, const double abseta) {
    int ptcorr_idx = -1;
    int abseta_idx = -1;
    for (int i = 0; i < ptcorr_bins.size()-1; i++) {
        if (ptcorr_bins[i] <= ptcorr && ptcorr < ptcorr_bins[i+1]) {
            ptcorr_idx = i;
            break;
        }
    }
    if (ptcorr_idx == -1) ptcorr_idx = ptcorr_bins.size()-2;

    for (int i = 0; i < abseta_bins.size()-1; i++) {
        if (abseta_bins[i] <= abseta && abseta < abseta_bins[i+1]) {
            abseta_idx = i;
            break;
        }
    }

    TString etaBin;
    if (abseta_idx == 0) etaBin = "EB1";
    else if (abseta_idx == 1) etaBin = "EB2";
    else if (abseta_idx == 2) etaBin = "EE";
    else etaBin = "EE";

    return TString::Format("ptcorr_%dto%d_%s",
                          static_cast<int>(ptcorr_bins[ptcorr_idx]),
                          static_cast<int>(ptcorr_bins[ptcorr_idx+1]),
                          etaBin.Data());
}

float MeasFakeRateV4::getJetPtCut(const TString& selection) {
    if (selection.Contains("MotherJetPt_Up"))
        return 60.0;
    else if (selection.Contains("MotherJetPt_Down"))
        return (leptonType == LeptonType::MUON) ? 20.0 : 30.0;
    else
        return 40.0;
}

template<typename T>
int MeasFakeRateV4::getMotherJetFlavour(const T& lep, const RVec<Jet>& allJets) {
    if (IsDATA) return -999;  // Not applicable for data

    short jetIdx = lep.JetIdx();
    if (jetIdx < 0) return -1;  // No associated jet

    for (const auto& jet : allJets) {
        if (jet.OriginalIndex() == jetIdx) {
            return (jet.genJetIdx() < 0) ? -1 : jet.hadronFlavour();
        }
    }
    return -1;  // Jet not found
}

// Explicit template instantiations
template int MeasFakeRateV4::getMotherJetFlavour<Muon>(const Muon&, const RVec<Jet>&);
template int MeasFakeRateV4::getMotherJetFlavour<Electron>(const Electron&, const RVec<Jet>&);

TString MeasFakeRateV4::getFlavourSubdir(int flavour) {
    if (flavour == 5) return "bjet";
    if (flavour == 4) return "cjet";
    if (flavour == 0) return "ljet";
    return "pujet";  // 0, -1, or any other value
}

void MeasFakeRateV4::fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst) {
    if (syst != "Central") return;
    TString channelStr = channelToString(channel);
    if (channelStr == "NONE") channelStr = "PreSel";

    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/%s/cutflow", channelStr.Data(), ID.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}
