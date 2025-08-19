#include "MeasFakeRate.h"

MeasFakeRate::MeasFakeRate() : leptonType(LeptonType::MUON), currentID("") {}

MeasFakeRate::~MeasFakeRate() {}

void MeasFakeRate::initializeAnalyzer() {
    // Set flags
    MeasFakeMu8 = HasFlag("MeasFakeMu8");
    MeasFakeMu17 = HasFlag("MeasFakeMu17");
    MeasFakeEl8 = HasFlag("MeasFakeEl8");
    MeasFakeEl12 = HasFlag("MeasFakeEl12");
    MeasFakeEl23 = HasFlag("MeasFakeEl23");
    MeasFakeMu = MeasFakeMu8 || MeasFakeMu17;
    MeasFakeEl = MeasFakeEl8 || MeasFakeEl12 || MeasFakeEl23;
    RunSyst = HasFlag("RunSyst");

    // Determine lepton type and binning
    if (MeasFakeMu) {
        leptonType = LeptonType::MUON;
        ptcorr_bins = {10., 15., 20., 30., 50., 100., 200.};
        abseta_bins = {0., 0.9, 1.6, 2.4};
    } else if (MeasFakeEl) {
        leptonType = LeptonType::ELECTRON;
        ptcorr_bins = {15., 20., 25., 35., 50., 100., 200.};
        abseta_bins = {0., 0.8, 1.479, 2.5};
    } else {
        throw std::runtime_error("[MeasFakeRate::initializeAnalyzer] No lepton type specified by flags");
    }

    // Set IDs
    MuonIDs = new IDContainer("HcToWATight", "HcToWALoose");
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
    }

    // Initialize correction
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    // Initialize SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/FakeRateSystematics.yaml", MCSample, DataEra);
    }

}

void MeasFakeRate::executeEvent() {
    Event ev = GetEvent();
    
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    
    RVec<Muon> rawMuons = GetAllMuons();
    RVec<Electron> rawElectrons = GetAllElectrons();
    
    if (!ev.PassTrigger(isoSglLepTrig)) return;
    
    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Loop over IDs (loose and tight)
    RVec<TString> IDs = {"loose", "tight"};
    
    for (const auto& ID : IDs) {
        currentID = ID;
        
        if (!IsDATA && RunSyst && systHelper) {
            // Process Central objects and weight-only systematics
            RecoObjects centralObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, centralObjects, ID, "Central");
            
            if (selectedChannel != Channel::NONE) {
                WeightInfo centralWeights = getWeights(selectedChannel, ID, ev, centralObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, centralObjects, centralWeights, "Central");
                
                // Process weight-only systematics using Central objects
                processWeightOnlySystematics(selectedChannel, ID, ev, centralObjects, genParts);
                
                // Process systematics requiring evtLoopAgain (object variations)
                for (const auto& syst : *systHelper) {
                    TString systName = syst.iter_name;
                    
                    // Skip Central and weight-only systematics
                    if (systName == "Central" || !systHelper->findSystematic(syst.syst_name)->evtLoopAgain) continue;
                    
                    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, systName);
                    Channel systChannel = selectEvent(ev, recoObjects, ID, systName);
                    
                    if (systChannel != Channel::NONE) {
                        WeightInfo weights = getWeights(systChannel, ID, ev, recoObjects, genParts, systName);
                        fillObjects(systChannel, ID, recoObjects, weights, systName);
                    }
                }
            }
        } else {
            // Process only Central for DATA or when systematics are off
            RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, recoObjects, ID, "Central");
            
            if (selectedChannel != Channel::NONE) {
                WeightInfo weights = getWeights(selectedChannel, ID, ev, recoObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, recoObjects, weights, "Central");
            }
        }
    }
}

MeasFakeRate::RecoObjects MeasFakeRate::defineObjects(Event& ev, 
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
    }
    
    // Get MET
    Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);
    METv = ApplyTypeICorrection(METv, allJets, allElectrons, allMuons);
    
    // Sort objects by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select objects based on ID
    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    
    RVec<Muon> looseMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Electron> looseElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 15., 2.5);
    
    RVec<Muon> tightMuons;
    RVec<Electron> tightElectrons;
    
    if (ID == "loose") {
        tightMuons = looseMuons;
        tightElectrons = looseElectrons;
    } else if (ID == "tight") {
        tightMuons = SelectMuons(allMuons, MuonIDs->GetID("tight"), 10., 2.4);
        tightElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("tight"), 15., 2.5);
    }
    
    // Jet selection with potential selection variations
    const float jetPtCut = getJetPtCut(syst);
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", jetPtCut, max_jeteta);
    if (Run == 2) tightJets = SelectJets(tightJets, "loosePuId", jetPtCut, max_jeteta);
    RVec<Jet> tightJets_vetoLep = JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4);
    
    // B-jet selection
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& jet : tightJets_vetoLep) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }
    
    RecoObjects objects;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.vetoMuons = vetoMuons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.vetoElectrons = vetoElectrons;
    objects.tightJets = tightJets;
    objects.tightJets_vetoLep = tightJets_vetoLep;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv = METv;
    
    return objects;
}

MeasFakeRate::Channel MeasFakeRate::selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst) {
    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;
    
    if (leptonType == LeptonType::MUON) {
        bool singleMu = (muons.size() == 1 && recoObjects.looseMuons.size() == 1 && 
                        electrons.size() == 0 && recoObjects.looseElectrons.size() == 0);
        bool doubleMu = (muons.size() == 2 && recoObjects.looseMuons.size() == 2 && 
                        electrons.size() == 0 && recoObjects.looseElectrons.size() == 0);
        
        if (!singleMu && !doubleMu) return Channel::NONE;
        if (muons[0].Pt() <= trigSafePtCut) return Channel::NONE;
        if (jets.size() == 0) return Channel::NONE;
        
        // RequireHeavyTag selection variation
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        
        if (singleMu) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(muons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            return Channel::INCLUSIVE;
        } else { // doubleMu
            Particle ZCand = muons[0] + muons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            return Channel::ZENRICHED;
        }
    } else { // ELECTRON
        bool singleEl = (electrons.size() == 1 && recoObjects.vetoElectrons.size() == 1 && 
                        muons.size() == 0 && recoObjects.vetoMuons.size() == 0);
        bool doubleEl = (electrons.size() == 2 && recoObjects.vetoElectrons.size() == 2 && 
                        muons.size() == 0 && recoObjects.vetoMuons.size() == 0);
        
        if (!singleEl && !doubleEl) return Channel::NONE;
        if (electrons[0].Pt() <= trigSafePtCut) return Channel::NONE;
        if (jets.size() == 0) return Channel::NONE;
        
        // RequireHeavyTag selection variation
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        
        if (singleEl) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(electrons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            return Channel::INCLUSIVE;
        } else { // doubleEl
            Particle ZCand = electrons[0] + electrons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            return Channel::ZENRICHED;
        }
    }
    
    return Channel::NONE;
}

MeasFakeRate::WeightInfo MeasFakeRate::getWeights(const Channel& channel,
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
    
    if (!IsDATA) {
        Event& ev = const_cast<Event&>(event);
        weights.genWeight = MCweight() * ev.GetTriggerLumi("Full");
        
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
            weights.btagSF = myCorr->GetBTaggingReweightMethod1a(recoObjects.tightJets_vetoLep, 
                                                               JetTagging::JetFlavTagger::DeepJet,
                                                               JetTagging::JetFlavTaggerWP::Medium,
                                                               JetTagging::JetTaggingSFMethod::mujets,
                                                               MyCorrection::variation::nom);
        }
        
        // Pileup reweighting
        if (syst.Contains("PileupReweight")) {
            weights.pileupWeight = myCorr->GetPUWeight(ev.nTrueInt(), var);
        } else {
            weights.pileupWeight = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
        }
    }
    
    return weights;
}

void MeasFakeRate::fillObjects(const Channel& channel,
                               const TString& ID,
                               const RecoObjects& recoObjects, 
                               const WeightInfo& weights, 
                               const TString& syst) {
    
    float totalWeight = weights.genWeight * weights.prefireWeight * weights.pileupWeight * 
                       weights.topPtWeight * weights.muonRecoSF * weights.eleRecoSF * 
                       weights.btagSF;
    
    TString prefix = channelToString(channel) + "/" + ID + "/" + syst;
    
    // Fill lepton histograms
    if (leptonType == LeptonType::MUON) {
        const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
        for (const auto& mu : muons) {
            TString binName = findBin(mu.Pt(), fabs(mu.Eta()));
            FillHist(prefix + "/muons/" + binName + "/pt", mu.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/muons/" + binName + "/eta", mu.Eta(), totalWeight, 50, -2.5, 2.5);
            FillHist(prefix + "/muons/" + binName + "/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
        }
    } else {
        const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
        for (const auto& ele : electrons) {
            TString binName = findBin(ele.Pt(), fabs(ele.Eta()));
            FillHist(prefix + "/electrons/" + binName + "/pt", ele.Pt(), totalWeight, 300, 0., 300.);
            FillHist(prefix + "/electrons/" + binName + "/eta", ele.Eta(), totalWeight, 50, -2.5, 2.5);
            FillHist(prefix + "/electrons/" + binName + "/phi", ele.Phi(), totalWeight, 64, -3.2, 3.2);
        }
    }
    
    // Fill jet histograms
    FillHist(prefix + "/jets/n", recoObjects.tightJets_vetoLep.size(), totalWeight, 20, 0., 20.);
    for (const auto& jet : recoObjects.tightJets_vetoLep) {
        FillHist(prefix + "/jets/pt", jet.Pt(), totalWeight, 300, 0., 300.);
        FillHist(prefix + "/jets/eta", jet.Eta(), totalWeight, 50, -2.5, 2.5);
        FillHist(prefix + "/jets/phi", jet.Phi(), totalWeight, 64, -3.2, 3.2);
    }
    
    // Fill b-jet histograms
    FillHist(prefix + "/bjets/n", recoObjects.bjets.size(), totalWeight, 10, 0., 10.);
    for (const auto& bjet : recoObjects.bjets) {
        FillHist(prefix + "/bjets/pt", bjet.Pt(), totalWeight, 300, 0., 300.);
        FillHist(prefix + "/bjets/eta", bjet.Eta(), totalWeight, 50, -2.5, 2.5);
        FillHist(prefix + "/bjets/phi", bjet.Phi(), totalWeight, 64, -3.2, 3.2);
    }
    
    // Fill MET histograms
    FillHist(prefix + "/met/pt", recoObjects.METv.Pt(), totalWeight, 300, 0., 300.);
    FillHist(prefix + "/met/phi", recoObjects.METv.Phi(), totalWeight, 64, -3.2, 3.2);
}

void MeasFakeRate::processWeightOnlySystematics(const Channel& channel, const TString& ID, const Event& event, const RecoObjects& recoObjects, const RVec<Gen>& genParts) {
    // Get weight-only systematics from SystematicHelper
    const std::vector<std::string>& allWeightOnlySystematics = systHelper->getWeightOnlySystematics();
    
    // Process each weight-only systematic (Up and Down variations)
    for (const std::string& systName : allWeightOnlySystematics) {
        // Up variation
        TString systNameUp = systName + "_Up";
        WeightInfo weightsUp = getWeights(channel, ID, event, recoObjects, genParts, systNameUp);
        fillObjects(channel, ID, recoObjects, weightsUp, systNameUp);
        
        // Down variation
        TString systNameDown = systName + "_Down";
        WeightInfo weightsDown = getWeights(channel, ID, event, recoObjects, genParts, systNameDown);
        fillObjects(channel, ID, recoObjects, weightsDown, systNameDown);
    }
}

TString MeasFakeRate::findBin(const double ptcorr, const double abseta) {
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

double MeasFakeRate::getJetPtCut(const TString& selection) {
    if (selection.Contains("MotherJetPtUp")) return 25.0;
    if (selection.Contains("MotherJetPtDown")) return 15.0;
    return 20.0; // default
}

