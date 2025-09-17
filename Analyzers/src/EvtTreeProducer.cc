#include "EvtTreeProducer.h"

EvtTreeProducer::EvtTreeProducer(){}
EvtTreeProducer::~EvtTreeProducer(){}

void EvtTreeProducer::WriteHist(){
    GetOutfile()->cd();
    newtree->Write();
}

void EvtTreeProducer::initializeAnalyzer() {
    Run1E2Mu = HasFlag("Run1E2Mu");
    Run3Mu = HasFlag("Run3Mu");

    MuonIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    ElectronIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    
    // initialize tree
    GetOutfile()->cd();
    newtree = new TTree("Events","Events");

    // muons
    newtree->Branch("nMuons", &nMuons);
    newtree->Branch("MuonPtColl", MuonPtColl, "MuonPtColl[nMuons]/F");
    newtree->Branch("MuonEtaColl", MuonEtaColl, "MuonEtaColl[nMuons]/F");
    newtree->Branch("MuonPhiColl", MuonPhiColl, "MuonPhiColl[nMuons]/F");
    newtree->Branch("MuonMassColl", MuonMassColl, "MuonMassColl[nMuons]/F");
    newtree->Branch("MuonChargeColl", MuonChargeColl, "MuonChargeColl[nMuons]/I");
    newtree->Branch("MuonLabelColl", MuonLabelColl, "MuonLabelColl[nMuons]/O");

    // electrons
    newtree->Branch("nElectrons", &nElectrons);
    newtree->Branch("ElectronPtColl", ElectronPtColl, "ElectronPtColl[nElectrons]/F");
    newtree->Branch("ElectronEtaColl", ElectronEtaColl, "ElectronEtaColl[nElectrons]/F");
    newtree->Branch("ElectronPhiColl", ElectronPhiColl, "ElectronPhiColl[nElectrons]/F");
    newtree->Branch("ElectronMassColl", ElectronMassColl, "ElectronMassColl[nElectrons]/F");
    newtree->Branch("ElectronChargeColl", ElectronChargeColl, "ElectronChargeColl[nElectrons]/I");
    newtree->Branch("ElectronLabelColl", ElectronLabelColl, "ElectronLabelColl[nElectrons]/O");

    // jets
    newtree->Branch("nJets", &nJets);
    newtree->Branch("JetPtColl", JetPtColl, "JetPtColl[nJets]/F");
    newtree->Branch("JetEtaColl", JetEtaColl, "JetEtaColl[nJets]/F");
    newtree->Branch("JetPhiColl", JetPhiColl, "JetPhiColl[nJets]/F");
    newtree->Branch("JetMassColl", JetMassColl, "JetMassColl[nJets]/F");
    newtree->Branch("JetChargeColl", JetChargeColl, "JetChargeColl[nJets]/I");
    newtree->Branch("JetBtagScoreColl", JetBtagScoreColl, "JetBtagScoreColl[nJets]/F");
    newtree->Branch("JetIsBtaggedColl", JetIsBtaggedColl, "JetIsBtaggedColl[nJets]/O");
    newtree->Branch("JetLabelColl", JetLabelColl, "JetLabelColl[nJets]/O");

    // METv
    newtree->Branch("METvPt", &METvPt);
    newtree->Branch("METvPhi", &METvPhi);

    // Correction
    if (IsDATA) {
        std::cerr << "EvtTreeProducer is not supported for data" << std::endl;
        exit(EXIT_FAILURE);
    }
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample ,IsDATA);
}

void EvtTreeProducer::executeEvent(){
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    RVec<Electron> rawElectrons = GetAllElectrons();
    
    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, "Central");
    Channel channel = selectEvent(ev, recoObjects, "Central");
    if (channel == Channel::NONE) return;

    nMuons = recoObjects.vetoMuons.size();
    for (std::size_t i = 0; i < nMuons; i++) {
        MuonPtColl[i] = recoObjects.vetoMuons[i].Pt();
        MuonEtaColl[i] = recoObjects.vetoMuons[i].Eta();
        MuonPhiColl[i] = recoObjects.vetoMuons[i].Phi();
        MuonMassColl[i] = recoObjects.vetoMuons[i].M();
        MuonChargeColl[i] = recoObjects.vetoMuons[i].Charge();
    }

    nElectrons = recoObjects.vetoElectrons.size();
    for (std::size_t i = 0; i < nElectrons; i++) {
        ElectronPtColl[i] = recoObjects.vetoElectrons[i].Pt();
        ElectronEtaColl[i] = recoObjects.vetoElectrons[i].Eta();
        ElectronPhiColl[i] = recoObjects.vetoElectrons[i].Phi();
        ElectronMassColl[i] = recoObjects.vetoElectrons[i].M();
    }

    nJets = recoObjects.jets.size();
    nBJets = recoObjects.bjets.size();
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (std::size_t i = 0; i < nJets; i++) {
        JetPtColl[i] = recoObjects.jets[i].Pt();
        JetEtaColl[i] = recoObjects.jets[i].Eta();
        JetPhiColl[i] = recoObjects.jets[i].Phi();
        JetMassColl[i] = recoObjects.jets[i].M();
        JetChargeColl[i] = recoObjects.jets[i].Charge();
        JetBtagScoreColl[i] = recoObjects.jets[i].GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        JetIsBtaggedColl[i] = recoObjects.jets[i].GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet) > wp;
    }

    METvPt = recoObjects.METv.Pt();
    METvPhi = recoObjects.METv.Phi();
    newtree->Fill();
}

EvtTreeProducer::RecoObjects EvtTreeProducer::defineObjects(const Event& ev, 
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
    
    // Sort objects in pt order
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    RVec<Muon> vetoMuons = SelectMuons(allMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> looseMuons = SelectMuons(vetoMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(looseMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<Electron> vetoElectrons = SelectElectrons(allElectrons, ElectronIDs->GetID("loose"), 10., 2.5);
    RVec<Electron> looseElectrons = SelectElectrons(vetoElectrons, ElectronIDs->GetID("loose"), 15., 2.5);
    RVec<Electron> tightElectrons = SelectElectrons(looseElectrons, ElectronIDs->GetID("tight"), 15., 2.5);
    
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> jets = SelectJets(allJets, "tight", 20., max_jeteta);
    if (Run == 2) {
        jets = SelectJets(jets, "loosePuId", 20., max_jeteta);

    }
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);

    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (auto& jet : jets) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
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

EvtTreeProducer::Channel EvtTreeProducer::selectEvent(const Event& ev, const RecoObjects& recoObjects, const TString& syst) {
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Muon>& looseMuons = recoObjects.looseMuons;
    const RVec<Muon>& tightMuons = recoObjects.tightMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Electron>& looseElectrons = recoObjects.looseElectrons;
    const RVec<Electron>& tightElectrons = recoObjects.tightElectrons;
    const RVec<Jet>& jets = recoObjects.jets;
    const RVec<Jet>& bjets = recoObjects.bjets;
    
    bool is1E2Mu = looseMuons.size() == 2 && looseElectrons.size() == 1 && vetoMuons.size() == 2 && vetoElectrons.size() == 1;
    bool is3Mu = looseMuons.size() == 3 && looseElectrons.size() == 0 && vetoMuons.size() == 3 && vetoElectrons.size() == 0;

    if (Run1E2Mu && (!is1E2Mu)) return Channel::NONE;
    if (Run3Mu && (!is3Mu)) return Channel::NONE;

    if (is1E2Mu) {
        // Reduced baseline for 1E2Mu channel
        // only require OS muon pair mass > 12 GeV
        // 60 < pair mass < 120 GeV
        // at least two jets 
        const Muon& mu1 = looseMuons[0];
        const Muon& mu2 = looseMuons[1];
        if (! (mu1.Charge() + mu2.Charge() == 0)) return Channel::NONE;
        Particle pair = mu1 + mu2;
        if (! (pair.M() > 12.)) return Channel::NONE;
        if (! (60 < pair.M() && pair.M() < 120)) return Channel::NONE;
        if (! (jets.size() >= 2)) return Channel::NONE;
        return Channel::SR1E2MU;
    } else { 
        // Reduced baseline for 3Mu channel
        // Exist OS pairs
        // all OS pairs mass > 12 GeV
        // at least one 60 < OS pair mass < 120 GeV
        // at least two jets
        const Muon& mu1 = looseMuons[0];
        const Muon& mu2 = looseMuons[1];
        const Muon& mu3 = looseMuons[2];
        if (! (mu1.Charge() + mu2.Charge() + mu3.Charge() == 1)) return Channel::NONE;

        // Make OS pairs
        Particle pair1, pair2;
        if (mu1.Charge() == mu2.Charge()) {
            pair1 = mu1 + mu3;
            pair2 = mu2 + mu3;
        } else if (mu1.Charge() == mu3.Charge()) {
            pair1 = mu1 + mu2;
            pair2 = mu2 + mu3;
        } else {
            pair1 = mu1 + mu2;
            pair2 = mu1 + mu3;
        }
        if (! (pair1.M() > 12.)) return Channel::NONE;
        if (! (pair2.M() > 12.)) return Channel::NONE;
        if (! ((60 < pair1.M() && pair1.M() < 120) || (60 < pair2.M() && pair2.M() < 120))) return Channel::NONE;
        if (! (jets.size() >= 2)) return Channel::NONE;
        return Channel::SR3MU;
    }
}
