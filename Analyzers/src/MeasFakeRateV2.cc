#include "MeasFakeRateV2.h"

MeasFakeRateV2::MeasFakeRateV2() : leptonType(LeptonType::NONE) {}

MeasFakeRateV2::~MeasFakeRateV2() {}

void MeasFakeRateV2::initializeAnalyzer() {
    // Set flags
    MeasFakeMu8 = HasFlag("MeasFakeMu8");
    MeasFakeMu17 = HasFlag("MeasFakeMu17");
    MeasFakeEl8 = HasFlag("MeasFakeEl8");
    MeasFakeEl12 = HasFlag("MeasFakeEl12");
    MeasFakeEl23 = HasFlag("MeasFakeEl23");
    RunSyst = HasFlag("RunSyst");

    // Determine lepton type and binning
    if (MeasFakeMu8 || MeasFakeMu17) {
        leptonType = LeptonType::MUON;
        ptcorr_bins = {10., 15., 20., 30., 50., 100., 200.};
        if (Run == 3) {
            eta_bins = {-2.4, -1.6, -0.9, 0., 0.9, 1.6, 2.4};
        } else {
            abseta_bins = {0., 0.9, 1.6, 2.4};
        }
    } else if (MeasFakeEl8 || MeasFakeEl12 || MeasFakeEl23) {
        leptonType = LeptonType::ELECTRON;
        ptcorr_bins = {10., 15., 20., 25., 35., 50., 100., 200.};
        if (Run == 3) {
            eta_bins = {-2.5, -1.479, -0.8, 0., 0.8, 1.479, 2.5};
        } else {
            abseta_bins = {0., 0.8, 1.479, 2.5};
        }
    } else {
        throw std::runtime_error("[MeasFakeRateV2::initializeAnalyzer] No lepton type specified by flags");
    }

    // Set IDs
    MuonIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    ElectronIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));

    // Set triggers
    if (MeasFakeMu8) {
        isoSglLepTrig = "HLT_Mu8_TrkIsoVVL";
        trigSafePtCut = 10.;
    } else if (MeasFakeMu17) {
        isoSglLepTrig = "HLT_Mu17_TrkIsoVVL";
        trigSafePtCut = 20.;
    } else if (MeasFakeEl8) {
        isoSglLepTrig = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30";
        trigSafePtCut = 10.;
    } else if (MeasFakeEl12) {
        isoSglLepTrig = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30";
        trigSafePtCut = 15.;
    } else if (MeasFakeEl23) {
        isoSglLepTrig = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30";
        trigSafePtCut = 25.;
    } else {
        throw std::runtime_error("[MeasFakeRateV2::initializeAnalyzer] No trigger specified by userflags");
    }

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

void MeasFakeRateV2::executeEvent() {
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
    if (!ev.PassTrigger(isoSglLepTrig)) return;
    fillCutflow(CutStage::Trigger, Channel::NONE, "event", initialWeight, "Central");
    
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

MeasFakeRateV2::RecoObjects MeasFakeRateV2::defineObjects(Event& ev, 
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
    
    // Get MET
    Particle METv = ApplyTypeICorrection(ev.GetMETVector(Event::MET_Type::PUPPI), allJets, allElectrons, allMuons);
    
    // Sort objects by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select objects based on ID
    RVec<Muon> looseMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Electron> looseElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    
    RVec<Muon> tightMuons;
    RVec<Electron> tightElectrons;
    
    if (ID == "loose") {
        tightMuons = looseMuons;
        tightElectrons = looseElectrons;
    } else if (ID == "tight") {
        tightMuons = SelectMuons(allMuons, MuonIDs->GetID("tight"), 10., 2.4);
        tightElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("tight"), 10., 2.5);
    }
    
    // Jet selection with potential selection variations
    const float jetPtCut = getJetPtCut(syst);
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", jetPtCut, max_jeteta);
    if (Run == 2) {
        RVec<Jet> tightJets_vetoMap;
        for (const auto &jet : SelectJets(tightJets, "loosePuId", jetPtCut, max_jeteta)) {
            if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.push_back(jet);
        }
        tightJets = tightJets_vetoMap;
    }
    tightJets = JetsVetoLeptonInside(tightJets, looseElectrons, looseMuons, 0.4);

    // B-jet selection
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& jet : tightJets) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }
    
    RecoObjects objects;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.tightJets = tightJets;
    //objects.tightJets_vetoLep = tightJets_vetoLep;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv = METv;
    
    return objects;
}

MeasFakeRateV2::Channel MeasFakeRateV2::selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst) {
    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
    const RVec<Muon>& vetoMuons = recoObjects.looseMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.looseElectrons;
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
        
        if (! (muons[0].Pt() > trigSafePtCut)) return Channel::NONE;
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
        
        if (! (electrons[0].Pt() > trigSafePtCut)) return Channel::NONE;
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

MeasFakeRateV2::WeightInfo MeasFakeRateV2::getWeights(const Channel& channel,
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
        if (Run == 2) {
            const RVec<Jet>& jets = recoObjects.tightJets;
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

void MeasFakeRateV2::fillObjects(const Channel& channel,
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
    
    TString prefix = channelToString(channel) + "/" + ID + "/" + syst;
    const RVec<Jet> &jets = recoObjects.tightJets;
    const RVec<Jet> &bjets = recoObjects.bjets;
    const Particle METv = recoObjects.METv;
    // Fill histograms based on channel and lepton type
    if (channel == Channel::INCLUSIVE) {
        if (leptonType == LeptonType::MUON) {
            const Muon& mu = (ID == "loose") ? recoObjects.looseMuons[0] : recoObjects.tightMuons[0];
            float ptcorr = mu.Pt()*(1.+max(0., mu.MiniPFRelIso()-0.1));
            float eta_value = (Run == 3) ? mu.Eta() : fabs(mu.Eta());
            float mT = TMath::Sqrt(2.*mu.Pt()*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
            float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
            TString binName = getBinPrefix(ptcorr, eta_value);
                
            // Fill inclusive channel histograms
            FillHist(prefix + "/muon/pt", mu.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
            FillHist(prefix + "/muon/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/muon/ptcorr", ptcorr, totalWeight, ptcorr_bins);
            if (Run == 3) {
                FillHist(prefix + "/muon/eta", mu.Eta(), totalWeight, eta_bins);
            } else {
                FillHist(prefix + "/muon/abseta", fabs(mu.Eta()), totalWeight, abseta_bins);
            }
            FillHist(prefix + "/MT", mT, totalWeight, 600, 0., 300.);
            FillHist(prefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
            FillHist(prefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
            FillHist(prefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
            FillHist(prefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                
            // Fill binned histograms
            TString binnedPrefix = binName + "/" + channelToString(channel) + "/" + ID + "/" + syst;
            FillHist(binnedPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
            if (Run == 3) {
                FillHist(binnedPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
            } else {
                FillHist(binnedPrefix + "/muon/abseta", fabs(mu.Eta()), totalWeight, 24, 0., 2.4);
            }
            FillHist(binnedPrefix + "/MT", mT, totalWeight, 500, 0., 500.);
            FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
            FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 600, 0., 300.);
                
            // Fill subchannel - determine QCDEnriched or WEnriched
            TString subchannel = "";
            if (mT < 25. && METv.Pt() < 25.) {
                subchannel = "QCDEnriched";
            } else if (mT > 60.) {
                subchannel = "WEnriched";
            } else {
                return;
            }
                
            TString subchannelPrefix = binName + "/" + subchannel + "/" + ID + "/" + syst;
            FillHist(subchannelPrefix + "/muon/pt", mu.Pt(), totalWeight, 200, 0., 200.);
            FillHist(subchannelPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
            FillHist(subchannelPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
            
            if (Run == 3) {
                FillHist(subchannelPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
            } else {
                FillHist(subchannelPrefix + "/muon/abseta", fabs(mu.Eta()), totalWeight, 24, 0., 2.4);
            }
            FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
            FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
        } else { // ELECTRON
            const Electron& el = (ID=="loose") ? recoObjects.looseElectrons[0]: recoObjects.tightElectrons[0];
            float ptcorr = el.Pt()*(1.+max(0., el.MiniPFRelIso()-0.1));
            float eta_value = (Run == 3) ? el.scEta() : fabs(el.scEta());
            float mT = TMath::Sqrt(2.*el.Pt()*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
            float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
            TString binName = getBinPrefix(ptcorr, eta_value);
                
            // Fill inclusive channel histograms
            FillHist(prefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
            FillHist(prefix + "/electron/phi", el.Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
            if (Run == 3) {
                FillHist(prefix + "/electron/eta", el.scEta(), totalWeight, eta_bins);
            } else {
                FillHist(prefix + "/electron/abseta", fabs(el.scEta()), totalWeight, abseta_bins);
            }
            FillHist(prefix + "/MT", mT, totalWeight, 600, 0., 300.);
            FillHist(prefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
            FillHist(prefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
            FillHist(prefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
            FillHist(prefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                
            // Fill binned histograms
            TString binnedPrefix = binName+"/"+channelToString(channel)+"/"+ID+"/"+syst;
            FillHist(binnedPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
            if (Run == 3) {
                FillHist(binnedPrefix + "/electron/eta", el.scEta(), totalWeight, eta_bins);
            } else {
                FillHist(binnedPrefix + "/electron/abseta", fabs(el.scEta()), totalWeight, abseta_bins);
            }
            FillHist(binnedPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
            FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
            FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                
            // Fill subchannel - determine QCDEnriched or WEnriched  
            TString subchannel = "";
            if (mT < 25. && METv.Pt() < 25.) {
                subchannel = "QCDEnriched";
            } else if (mT > 60.) {
                subchannel = "WEnriched";
            } else {
                return;
            }
                
            TString subchannelPrefix = binName+"/"+subchannel+"/"+ID+"/"+syst;
            FillHist(subchannelPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
            FillHist(subchannelPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
            FillHist(subchannelPrefix + "/electron/ptcorr", ptcorr, totalWeight, 300, 0., 300.);
            if (Run == 3) {
                FillHist(subchannelPrefix + "/electron/eta", el.scEta(), totalWeight, 50, -2.5, 2.5);
            } else {
                FillHist(subchannelPrefix + "/electron/abseta", fabs(el.scEta()), totalWeight, 25, 0., 2.5);
            }
            FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
            FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
        }
    } else if (channel == Channel::ZENRICHED) {
        if (leptonType == LeptonType::MUON) {
            const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
            const Particle ZCand = muons[0] + muons[1];
            FillHist(prefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
            FillHist(prefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
            FillHist(prefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/muons/1/pt", muons[0].Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/muons/1/eta", muons[0].Eta(), totalWeight, 48, -2.4, 2.4);
            FillHist(prefix + "/muons/1/phi", muons[0].Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/muons/2/pt", muons[1].Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/muons/2/eta", muons[1].Eta(), totalWeight, 48, -2.4, 2.4);
            FillHist(prefix + "/muons/2/phi", muons[1].Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
            FillHist(prefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
        } else { // ELECTRON
            const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
            const Particle ZCand = electrons[0] + electrons[1];
            FillHist(prefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
            FillHist(prefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
            FillHist(prefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/electrons/1/pt", electrons[0].Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/electrons/1/scEta", electrons[0].scEta(), totalWeight, 50, -2.5, 2.5);
            FillHist(prefix + "/electrons/1/phi", electrons[0].Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/electrons/2/pt", electrons[1].Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/electrons/2/scEta", electrons[1].scEta(), totalWeight, 50, -2.5, 2.5);
            FillHist(prefix + "/electrons/2/phi", electrons[1].Phi(), totalWeight, 64, -3.2, 3.2);
            FillHist(prefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
            FillHist(prefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
        }
    }
}

TString MeasFakeRateV2::getBinPrefix(const double ptcorr, const double eta_value) {
    int ptcorr_idx = -1;
    int eta_idx = -1;
    
    for (int i = 0; i < ptcorr_bins.size()-1; i++) {
        if (ptcorr_bins[i] <= ptcorr && ptcorr < ptcorr_bins[i+1]) {
            ptcorr_idx = i;
            break;
        }
    }
    if (ptcorr_idx == -1) ptcorr_idx = ptcorr_bins.size()-2;
    
    if (Run == 3) {
        // Use signed eta for Run 3
        for (int i = 0; i < eta_bins.size()-1; i++) {
            if (eta_bins[i] <= eta_value && eta_value < eta_bins[i+1]) {
                eta_idx = i;
                break;
            }
        }
        if (eta_idx == -1) eta_idx = eta_bins.size()-2;
    } else {
        // Use absolute eta for Run 2
        double abseta = fabs(eta_value);
        for (int i = 0; i < abseta_bins.size()-1; i++) {
            if (abseta_bins[i] <= abseta && abseta < abseta_bins[i+1]) {
                eta_idx = i;
                break;
            }
        }
        if (eta_idx == -1) eta_idx = abseta_bins.size()-2;
    }
    
    TString etaBin;
    if (Run == 3) {
        // For Run 3 with signed eta, use different naming convention
        if (eta_idx == 0) etaBin = "EEm";      // -2.4/-2.5 to -1.6/-1.479
        else if (eta_idx == 1) etaBin = "EB2m"; // -1.6/-1.479 to -0.9/-0.8
        else if (eta_idx == 2) etaBin = "EB1m"; // -0.9/-0.8 to 0
        else if (eta_idx == 3) etaBin = "EB1p"; // 0 to 0.9/0.8
        else if (eta_idx == 4) etaBin = "EB2p"; // 0.9/0.8 to 1.6/1.479
        else if (eta_idx == 5) etaBin = "EEp";  // 1.6/1.479 to 2.4/2.5
        else etaBin = "EEp";
    } else {
        // For Run 2 with absolute eta
        if (eta_idx == 0) etaBin = "EB1";
        else if (eta_idx == 1) etaBin = "EB2";
        else if (eta_idx == 2) etaBin = "EE";
        else etaBin = "EE";
    }
    
    return TString::Format("ptcorr_%dto%d_%s", 
                          static_cast<int>(ptcorr_bins[ptcorr_idx]), 
                          static_cast<int>(ptcorr_bins[ptcorr_idx+1]), 
                          etaBin.Data());
}

float MeasFakeRateV2::getJetPtCut(const TString& selection) {
    if (selection.Contains("MotherJetPt_Up")) 
        return 60.0;
    else if (selection.Contains("MotherJetPt_Down")) 
        return (leptonType == LeptonType::MUON) ? 20.0 : 30.0;
    else 
        return 40.0;
}

void MeasFakeRateV2::fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst) {
    if (syst != "Central") return;
    TString channelStr = channelToString(channel);
    if (channelStr == "NONE") channelStr = "PreSel";
    
    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/%s/cutflow", channelStr.Data(), ID.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}

