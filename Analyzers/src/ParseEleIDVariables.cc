#include "ParseEleIDVariables.h"

ParseEleIDVariables::ParseEleIDVariables() {}
ParseEleIDVariables::~ParseEleIDVariables() {}

void ParseEleIDVariables::initializeAnalyzer() {
    Events = new TTree("Events", "Events");
    Events->Branch("genWeight", &genWeight);
    Events->Branch("nElectrons", &nElectrons);
    Events->Branch("pt", pt, "pt[nElectrons]/F");
    Events->Branch("scEta", scEta, "scEta[nElectrons]/F");
    Events->Branch("sieie", sieie, "sieie[nElectrons]/F");
    Events->Branch("lepType", lepType, "lepType[nElectrons]/I");
    Events->Branch("deltaEtaInSC", deltaEtaInSC, "deltaEtaInSC[nElectrons]/F");
    Events->Branch("deltaPhiInSeed", deltaPhiInSeed, "deltaPhiInSC[nElectrons]/F");
    Events->Branch("hoe", hoe, "hoe[nElectrons]/F");
    Events->Branch("ecalPFClusterIso", ecalPFClusterIso, "ecalPFClusterIso[nElectrons]/F");
    Events->Branch("hcalPFClusterIso", hcalPFClusterIso, "hcalPFClusterIso[nElectrons]/F");
    Events->Branch("rho", rho, "rho[nElectrons]/F");
    Events->Branch("dr03TkSumPt", dr03TkSumPt, "dr03TkSumPt[nElectrons]/F");
    Events->Branch("isMVANoIsoWP90", isMVANoIsoWP90, "isMVANoIsoWP90[nElectrons]/O");
    Events->Branch("isPOGMedium", isPOGMedium, "isPOGMedium[nElectrons]/O");
    Events->Branch("isPOGTight", isPOGTight, "isPOGTight[nElectrons]/O");
    Events->Branch("convVeto",  convVeto, "convVeto[nElectrons]/O");
    Events->Branch("lostHits", lostHits, "lostHits[nElectrons]/I");
    Events->Branch("dZ", dZ, "dZ[nElectrons]/F");
    Events->Branch("sip3d", sip3d, "sip3d[nElectrons]/F");
    Events->Branch("miniPFRelIso", miniPFRelIso, "miniPFRelIso[nElectrons]/F");
    Events->Branch("mvaNoIso", mvaNoIso, "mvaNoIso[nElectrons]/F");
    Events->Branch("nearestJetFlavour", nearestJetFlavour, "nearestJetFlavour[nElectrons]/I");
    Events->Branch("isEMuTrigMatched", isEMuTrigMatched, "isEMuTrigMatched[nElectrons]/O");
    Events->Branch("isIsoElTrigMatched", isIsoElTrigMatched, "isIsoElTrigMatched[nElectrons]/O");

    if (DataEra.Contains("2016")) {
        SglMuTriggers = {"HLT_IsoMu24", "HLT_IsoTkMu24"};
    } else if (DataEra == "2017") {
        SglMuTriggers = {"HLT_IsoMu27"};
    } else {
        SglMuTriggers = {"HLT_IsoMu24"};
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, MCSample, IsDATA);
}

void ParseEleIDVariables::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> jets = GetAllJets();
    if (!PassNoiseFilter(jets, ev)) return;

    RVec<Electron> electrons = GetElectrons("", 15., 2.5);
    RVec<Muon> muons = GetMuons("POGTight", 25., 2.4);
    RVec<Gen> truth = GetAllGens();
    RVec<TrigObj> trigObjs = GetAllTrigObjs();

    // Require event to pass EMu trigger
    // and hard muon for the tag
    if (! ev.PassTrigger(SglMuTriggers)) return;
    if (! (muons.size() == 1)) return;
    const auto &mu = muons.at(0);
    const float safePtCut = (DataEra.Contains("2017") ? 30. : 27.);
    if (! (mu.Pt() > safePtCut)) return;
    if (! PassSLT(mu, trigObjs)) return;
    if (! (electrons.size() > 0)) return;
    
    // Update branches
    genWeight = MCweight()*ev.GetTriggerLumi("Full")*GetL1PrefireWeight()*myCorr->GetPUWeight(ev.nTrueInt());
    nElectrons = electrons.size();
    for (int i = 0; i < nElectrons; i++) {
        const auto &el = electrons.at(i);
        pt[i] = el.Pt();
        scEta[i] = el.scEta();
        sieie[i] = el.sieie();
        deltaEtaInSC[i] = el.deltaEtaInSC();
        deltaPhiInSeed[i] = el.deltaPhiInSeed();
        hoe[i] = el.hoe();
        ecalPFClusterIso[i] = el.ecalPFClusterIso();
        hcalPFClusterIso[i] = el.hcalPFClusterIso();
        rho[i] = el.rho();
        dr03TkSumPt[i] = el.dr03TkSumPt();
        isMVANoIsoWP90[i] = el.PassID("POGMVANoIsoWP90");
        isPOGMedium[i] = el.PassID("POGMedium");
        isPOGTight[i] = el.PassID("POGTight");
        convVeto[i] = el.ConvVeto();
        lostHits[i] = el.LostHits();
        sip3d[i] = el.SIP3D();
        dZ[i] = el.dZ();
        miniPFRelIso[i] = el.MiniPFRelIso();
        mvaNoIso[i] = el.MvaNoIso();
        lepType[i] = GetLeptonType(el, truth);

        // Use jetIdx for efficient jet matching
        nearestJetFlavour[i] = -1; // Default: no jet match
        
        short jetIdx = el.JetIdx();
        FillHist("electronJetIdx", jetIdx, 1.0, 100, -10., 90.);
        
        if (jetIdx >= 0) {
            // Find the jet with matching original index
            const Jet* matchedJet = nullptr;
            for (const auto &jet : jets) {
                if (jet.OriginalIndex() == jetIdx) {
                    matchedJet = &jet;
                    break;
                }
            }
            if (matchedJet) {
                nearestJetFlavour[i] = matchedJet->genJetIdx() < 0 ? -1 : matchedJet->hadronFlavour();
            }
        }
        
        // Check trigger object matching for CaloIdL_TrackIdL_IsoVL filter
        isEMuTrigMatched[i] = PassEMT(el, trigObjs);
        isIsoElTrigMatched[i] = PassIsoElT(el, trigObjs);
    }
    Events->Fill();
}

void ParseEleIDVariables::WriteHist() {
    TFile* outfile = GetOutfile();
    Events->Write();
    outfile->Close();
}

bool ParseEleIDVariables::PassSLT(const Muon &mu, const RVec<TrigObj> &trigObjs) {
    const float trig_pt_cut = (DataYear == 2017 ? 27: 24);
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(3))) continue;
        if (! (trigObj.Pt() > trig_pt_cut)) continue;
        return true;
    }
    return false;
}

bool ParseEleIDVariables::PassIsoElT(const Electron &el, const RVec<TrigObj> &trigObjs) {
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isElectron()) continue;
        if (! (trigObj.DeltaR(el) < 0.3)) continue;
        if (! trigObj.hasBit(0)) continue;
        return true;
    }
    return false;
}

bool ParseEleIDVariables::PassEMT(const Electron &el, const RVec<TrigObj> &trigObjs){
    const int trigBit = (Run == 2) ? 5 : 6; // 5 for Run2, 6 for Run3
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isElectron()) continue;
        if (! (trigObj.DeltaR(el) < 0.3)) continue;
        if (! (trigObj.hasBit(trigBit))) continue;
        //if (! (trigObj.Pt() > pt_cut)) continue;
        return true;
    }
    return false;
}
