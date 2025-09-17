#include "DiLepton.h"

DiLepton::DiLepton() : channel(Channel::NONE), run_syst(false) {}

DiLepton::~DiLepton() {}

void DiLepton::initializeAnalyzer() {
    DiLeptonBase::initializeAnalyzer();

    // Determine channel
    if (RunDiMu) channel = Channel::DIMU;
    if (RunEMu) channel = Channel::EMU;
    if (RunSyst) run_syst = true;

    // Initialize SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/DiLeptonSystematics.yaml", MCSample, DataEra);
    }
}

void DiLepton::executeEvent() {
    Event ev = GetEvent();
    
    // Initial cutflow entry
    float initialWeight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    fillCutflow(CutStage::Initial, Channel::NONE, initialWeight, "Central");
    
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    fillCutflow(CutStage::NoiseFilter, Channel::NONE, initialWeight, "Central");

    RVec<Muon> rawMuons = GetAllMuons();
    if (!(RunNoVetoMap || PassVetoMap(rawJets, rawMuons, "jetvetomap"))) return; // For Run3, reject events if any jet is within the veto map
    fillCutflow(CutStage::VetoMap, Channel::NONE, initialWeight, "Central");

    RVec<Electron> rawElectrons = GetAllElectrons();
    
    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Use SystematicHelper for systematic variations
    if (!IsDATA && run_syst && systHelper) {
        // Step 1: Process Central objects and weight-only systematics
        RecoObjects centralObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, "Central");
        fillCutflow(CutStage::LeptonSelection, Channel::NONE, initialWeight, "Central");
        
        Channel selectedChannel = selectEvent(ev, centralObjects, "Central");
        
        if (selectedChannel != Channel::NONE) {
            fillCutflow(CutStage::Final, selectedChannel, initialWeight, "Central");
            // Fill Central with nominal weights
            WeightInfo centralWeights = getWeights(selectedChannel, ev, centralObjects, genParts, "Central");
            fillObjects(selectedChannel, centralObjects, centralWeights, "Central");

            // Process weight-only systematics using Central objects
            processWeightOnlySystematics(selectedChannel, ev, centralObjects, genParts);
        }
            
        // Process systematics requiring evtLoopAgain
        for (const auto& syst : *systHelper) {
            TString systName = syst.iter_name;
                
            // Skip Central (already processed) and weight-only systematics
            if (systName == "Central" || !systHelper->findSystematic(syst.syst_name)->evtLoopAgain) continue;
                
            // Define objects with systematic variation
            RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, systName);
            Channel systChannel = selectEvent(ev, recoObjects, systName);
                
            if (systChannel != Channel::NONE) {
                WeightInfo weights = getWeights(systChannel, ev, recoObjects, genParts, systName);
                fillObjects(systChannel, recoObjects, weights, systName);
            }
        }
    } else {
        // Process only Central for DATA or when systematics are off
        RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, "Central");
        fillCutflow(CutStage::LeptonSelection, Channel::NONE, initialWeight, "Central");
        
        Channel selectedChannel = selectEvent(ev, recoObjects, "Central");
        
        if (selectedChannel != Channel::NONE) {
            fillCutflow(CutStage::Final, selectedChannel, initialWeight, "Central");
            WeightInfo weights = getWeights(selectedChannel, ev, recoObjects, genParts, "Central");
            fillObjects(selectedChannel, recoObjects, weights, "Central");
        }
    }
}

DiLepton::RecoObjects DiLepton::defineObjects(Event& ev, 
                                             const RVec<Muon>& rawMuons, 
                                             const RVec<Electron>& rawElectrons, 
                                             const RVec<Jet>& rawJets,
                                             const RVec<GenJet>& genJets, 
                                             const TString& syst) {

    // Create copies of the raw objects
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;
    
    // Determine systematic variation type and direction
    MyCorrection::variation systVar = MyCorrection::variation::nom;
    TString systSource = "";
    
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
        systVar = syst.Contains("Up") ? MyCorrection::variation::up : MyCorrection::variation::down;
        systSource = "total"; // Use default JES source for simple JetEn systematics
        allJets = ScaleJets(allJets, variation, systSource);
    } else if (syst.Contains("JetRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allJets = SmearJets(allJets, genJets, variation);
    } else {
        // No scale variation
    }
    
    // Get MET from event and re-apply Type-I correction and XY correction
    Particle METv_default;
    if (syst.Contains("UnclusteredEn")) {
        Event::MET_Syst variation = syst.Contains("Up") ? Event::MET_Syst::UE_UP : Event::MET_Syst::UE_DOWN;
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI, variation);
    } else {
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    }
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);
    // XY correction is not recommended for PUPPI MET
    
    // Example for using CHS MET
    // UE variation should be done inside of ApplyTypeICorrection
    // Particle METv_default = ev.GetMETVector(Event::MET_Type::CHS);
    // Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);
    // Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons, MyCorrection::variation::up); // unclustered up
    // Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons, MyCorrection::variation::down); // unclustered down
    // METv = myCorr->METXYCorrection(METv, ev.run(), ev.nPV(), MyCorrection::XYCorrection_MetType::Type1CHSMET);

    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("tight"), 15., 2.5);
    
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", 20., max_jeteta);
    if (Run == 2) {
        tightJets = SelectJets(tightJets, "loosePuId", 20., max_jeteta);
        RVec<Jet> tightJets_vetoMap;
        for (const auto &jet: tightJets_vetoMap)
            if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.emplace_back(jet);
        if (!RunNoVetoMap) tightJets = tightJets_vetoMap;
    }
    RVec<Jet> tightJets_vetoLep = JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4);

    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (auto& jet : tightJets_vetoLep) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    RecoObjects objects;
    objects.vetoMuons = vetoMuons;
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.tightElectrons = tightElectrons;
    objects.tightJets = tightJets;
    objects.tightJets_vetoLep = tightJets_vetoLep;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv_default = METv_default;
    objects.METv = METv;

    return objects;
}

DiLepton::Channel DiLepton::selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& syst) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;

    float weight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    bool isDiMu = (tightMuons.size() == 2 && vetoMuons.size() == 2 && 
                   tightElectrons.size() == 0 && vetoElectrons.size() == 0);
    bool isEMu = (tightMuons.size() == 1 && vetoMuons.size() == 1 && 
                  tightElectrons.size() == 1 && vetoElectrons.size() == 1);

    if (channel == Channel::DIMU) {
        if (!isDiMu) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::DIMU, weight, syst);
    }
    if (channel == Channel::EMU) {
        if (!isEMu) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::EMU, weight, syst);
    }

    // DiMu selection
    if (channel == Channel::DIMU) {
        if (!ev.PassTrigger(DblMuTriggers)) return Channel::NONE;
        fillCutflow(CutStage::Trigger, Channel::DIMU, weight, syst);
        
        const Muon& mu1 = tightMuons[0];
        const Muon& mu2 = tightMuons[1];
        if (mu1.Pt() <= 20.) return Channel::NONE;
        if (mu2.Pt() <= 10.) return Channel::NONE;
        if (mu1.Charge() + mu2.Charge() != 0) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (pair.M() <= 50.) return Channel::NONE;
        fillCutflow(CutStage::KinematicCuts, Channel::DIMU, weight, syst);
        return Channel::DIMU;
    }
    // EMu selection
    else if (channel == Channel::EMU) {
        if (!ev.PassTrigger(EMuTriggers)) return Channel::NONE;
        fillCutflow(CutStage::Trigger, Channel::EMU, weight, syst);
        
        const Muon& mu = tightMuons[0];
        const Electron& ele = tightElectrons[0];
        if (!((mu.Pt() > 25. && ele.Pt() > 15.) || (mu.Pt() > 10. && ele.Pt() > 25.))) return Channel::NONE;
        if (mu.Charge() + ele.Charge() != 0) return Channel::NONE;
        if (mu.DeltaR(ele) <= 0.4) return Channel::NONE;
        fillCutflow(CutStage::KinematicCuts, Channel::EMU, weight, syst);
        
        if (jets.size() < 2) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::EMU, weight, syst);
        
        if (bjets.size() < 1) return Channel::NONE;
        fillCutflow(CutStage::BjetRequirements, Channel::EMU, weight, syst);
        
        return Channel::EMU;
    }

    return Channel::NONE;
}

void DiLepton::processWeightOnlySystematics(const DiLepton::Channel& channel,
                                           const Event& event,
                                           const RecoObjects& recoObjects, 
                                           const RVec<Gen>& genParts) {
    // Get weight-only systematics from SystematicHelper and filter for channel
    const std::vector<std::string>& allWeightOnlySystematics = systHelper->getWeightOnlySystematics();
    std::vector<std::string> channelFilteredSystematics;
    
    for (const std::string& systName : allWeightOnlySystematics) {
        // Filter channel-specific systematics
        bool includeSystematic = true;
        if (systName == "DblMuTrigSF" && channel != Channel::DIMU) includeSystematic = false;
        if (systName == "EMuTrigSF" && channel != Channel::EMU) includeSystematic = false;
        if (systName == "ElectronIDSF" && channel != Channel::EMU) includeSystematic = false;
        
        if (includeSystematic) {
            channelFilteredSystematics.push_back(systName);
        }
    }

    // Process each weight-only systematic (Up and Down variations)
    for (const std::string& systName : channelFilteredSystematics) {
        // Up variation
        TString systNameUp = systName + "_Up";
        WeightInfo weightsUp = getWeights(channel, event, recoObjects, genParts, systNameUp);
        fillObjects(channel, recoObjects, weightsUp, systNameUp);
        
        // Down variation
        TString systNameDown = systName + "_Down";
        WeightInfo weightsDown = getWeights(channel, event, recoObjects, genParts, systNameDown);
        fillObjects(channel, recoObjects, weightsDown, systNameDown);

        // Check not implemented case
        TString systNameNotImplemented = systName + "_NotImplemented";
        WeightInfo centralWeights = getWeights(channel, event, recoObjects, genParts, "Central");
        fillObjects(channel, recoObjects, centralWeights, systNameNotImplemented);
    }
}

DiLepton::WeightInfo DiLepton::getWeights(const DiLepton::Channel& channel,
                                         const Event& event,
                                         const RecoObjects& recoObjects, 
                                         const RVec<Gen>& genParts, 
                                         const TString& syst) {
    WeightInfo weights;
    
    if (IsDATA) {
        weights.genWeight = 1.;
        weights.prefireWeight = 1.;
        weights.pileupWeight = 1.;
        weights.topPtWeight = 1.;
        weights.muonRecoSF = 1.;
        weights.muonIDSF = 1.;
        weights.eleRecoSF = 1.;
        weights.eleIDSF = 1.;
        weights.trigSF = 1.;
        weights.pileupIDSF = 1.;
        weights.btagSF = 1.;
        return weights;
    }

    // SystematicHelper handles systematic validation
    weights.genWeight = MCweight() * event.GetTriggerLumi("Full");

    // Use SystematicHelper for weight variations or nominal values
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

    // Top pT reweighting
    weights.topPtWeight = 1.;
    if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ")) {
        RVec<Gen> genParts = GetAllGens();
        weights.topPtWeight = myCorr->GetTopPtReweight(genParts);
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
    if (channel == Channel::EMU) {
        if (syst.Contains("EMuTrigSF")) {
            weights.trigSF = myCorr->GetEMuTriggerSF(electrons, muons, var);
        } else {
            weights.trigSF = myCorr->GetEMuTriggerSF(electrons, muons);
        }
    } else if (channel == Channel::DIMU) {
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
        const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
        const RVec<GenJet>& genJets = recoObjects.genJets;
        unordered_map<int, int> matched_idx = GenJetMatching(jets, genJets, fixedGridRhoFastjetAll, 0.4, 10.);
        if (syst.Contains("PileupJetIDSF")) {
            weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", var);
        } else {
            weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", MyCorrection::variation::nom);
        }
    } else {
        weights.pileupIDSF = 1.;
    }

    // B-tagging scale factor (EMu channel only)
    if (channel == Channel::EMU) {
        const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
        JetTagging::JetFlavTagger tagger = JetTagging::JetFlavTagger::DeepJet;
        JetTagging::JetFlavTaggerWP wp = JetTagging::JetFlavTaggerWP::Medium;
        JetTagging::JetTaggingSFMethod method = JetTagging::JetTaggingSFMethod::mujets;
        //weights.btagSF = myCorr->GetBTaggingSF(jets, tagger, wp, method);
        string source="central";
        if (syst.Contains("HFcorr")) source = "hf_corr";
        if (syst.Contains("HFuncorr")) source = "hf_uncorr";
        if (syst.Contains("LFcorr")) source = "lf_corr";
        if (syst.Contains("LFuncorr")) source = "lf_uncorr";
        weights.btagSF = myCorr->GetBTaggingReweightMethod1a(jets, tagger, wp, method, var, source);
    } else {
        weights.btagSF = 1.;
    }

    return weights;
}

void DiLepton::fillObjects(const DiLepton::Channel& channel, const RecoObjects& recoObjects, 
                          const WeightInfo& weights, const TString& syst) {
    const RVec<Muon>& muons = recoObjects.tightMuons;
    const RVec<Electron>& electrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets_vetoLep;
    const RVec<Jet>& bjets = recoObjects.bjets;
    const Particle& METv_default = recoObjects.METv_default;
    const Particle& METv = recoObjects.METv;

    TString channelStr = channelToString(channel);
    float weight = 1.;

    if (!IsDATA) {
        weight = weights.genWeight * weights.prefireWeight * weights.pileupWeight * 
                weights.topPtWeight * weights.muonRecoSF * weights.muonIDSF * 
                weights.eleRecoSF * weights.eleIDSF * weights.trigSF * 
                weights.pileupIDSF * weights.btagSF;

        if (syst.Contains("NotImplemented")) {
            if (syst.Contains("L1Prefire")) weight /= weights.prefireWeight;
            if (syst.Contains("PileupReweight")) weight /= weights.pileupWeight;
            if (syst.Contains("TopPtReweight")) weight /= weights.topPtWeight;
            if (syst.Contains("PileupJetIDSF")) weight /= weights.pileupIDSF;
            if (syst.Contains("MuonIDSF")) weight /= (weights.muonIDSF*weights.trigSF);
            if (syst.Contains("ElectronIDSF")) weight /= (weights.eleIDSF*weights.trigSF);
            if (syst.Contains("DblMuTrigSF")) weight /= weights.trigSF;
            if (syst.Contains("EMuTrigSF")) weight /= weights.trigSF;
            if (syst.Contains("BtagSF")) weight /= weights.btagSF;
        }


        FillHist(Form("%s/%s/weights/genWeight", channelStr.Data(), syst.Data()), weights.genWeight, 1., 200, -10000, 10000.);
        FillHist(Form("%s/%s/weights/prefireWeight", channelStr.Data(), syst.Data()), weights.prefireWeight, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/pileupWeight", channelStr.Data(), syst.Data()), weights.pileupWeight, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/topPtWeight", channelStr.Data(), syst.Data()), weights.topPtWeight, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/muonRecoSF", channelStr.Data(), syst.Data()), weights.muonRecoSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/muonIDSF", channelStr.Data(), syst.Data()), weights.muonIDSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/eleRecoSF", channelStr.Data(), syst.Data()), weights.eleRecoSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/eleIDSF", channelStr.Data(), syst.Data()), weights.eleIDSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/trigSF", channelStr.Data(), syst.Data()), weights.trigSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/pileupIDSF", channelStr.Data(), syst.Data()), weights.pileupIDSF, 1., 100, -5., 5.);
        FillHist(Form("%s/%s/weights/btagSF", channelStr.Data(), syst.Data()), weights.btagSF, 1., 100, -5., 5.);
    }

    // Fill muon histograms
    for (size_t idx = 0; idx < muons.size(); ++idx) {
        const Muon& mu = muons.at(idx);
        FillHist(Form("%s/%s/muons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), mu.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/muons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(Form("%s/%s/muons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/muons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), mu.M(), weight, 10, 0., 1.);
    }

    // Fill electron histograms
    for (size_t idx = 0; idx < electrons.size(); ++idx) {
        const Electron& ele = electrons.at(idx);
        FillHist(Form("%s/%s/electrons/%zu/pt", channelStr.Data(), syst.Data(), idx+1), ele.Pt(), weight, 300, 0., 300.);
        FillHist(Form("%s/%s/electrons/%zu/eta", channelStr.Data(), syst.Data(), idx+1), ele.Eta(), weight, 50, -2.5, 2.5);
        FillHist(Form("%s/%s/electrons/%zu/phi", channelStr.Data(), syst.Data(), idx+1), ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(Form("%s/%s/electrons/%zu/mass", channelStr.Data(), syst.Data(), idx+1), ele.M(), weight, 100, 0., 1.);
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
}

void DiLepton::fillCutflow(CutStage stage, const Channel& channel, float weight, const TString& syst) {
    if (syst != "Central") return;
    TString channelStr = channelToString(channel);
    if (channelStr == "NONE") channelStr = "ALL";
    
    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/cutflow", channelStr.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}
