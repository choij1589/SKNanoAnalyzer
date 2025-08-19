#include "ParseMuIDVariables.h"

ParseMuIDVariables::ParseMuIDVariables() {}
ParseMuIDVariables::~ParseMuIDVariables() {}

void ParseMuIDVariables::initializeAnalyzer() {
    Events = new TTree("Events", "Events");
    Events->Branch("genWeight", &genWeight);
    Events->Branch("nMuons", &nMuons);
    Events->Branch("pt", pt, "pt[nMuons]/F");
    Events->Branch("eta", eta, "eta[nMuons]/F");
    Events->Branch("lepType", lepType, "lepType[nMuons]/I");
    Events->Branch("nearestJetFlavour", nearestJetFlavour, "nearestJetFlavour[nMuons]/I");
    Events->Branch("isPOGMediumId", isPOGMediumId, "isPOGMediumId[nMuons]/F");
    Events->Branch("dZ", dZ, "dZ[nMuons]/F");
    Events->Branch("sip3d", sip3d, "sip3d[nMuons]/F");
    Events->Branch("tkRelIso", tkRelIso, "tkRelIso[nMuons]/F");
    Events->Branch("miniPFRelIso", miniPFRelIso, "miniPFRelIso[nMuons]/F");
    Events->Branch("isEMuTrigMatched", isEMuTrigMatched, "isEMuTrigMatched[nMuons]/O");
    Events->Branch("isIsoMuTrigMatched", isIsoMuTrigMatched, "isIsoMuTrigMatched[nMuons]/O");

    if (DataEra == "2016preVFP") {
        SglElTriggers = {"HLT_Ele27_WPTight_Gsf"};
    } else if (DataEra == "2016postVFP") {
        SglElTriggers = {"HLT_Ele27_WPTight_Gsf"};
    } else if (DataEra == "2017") {
        SglElTriggers = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG"};
    } else if (DataEra == "2018") {
        SglElTriggers = {"HLT_Ele32_WPTight_Gsf"};
    } else if (Run == 3) {
        SglElTriggers = {"HLT_Ele30_WPTight_Gsf"};
    } else {
        cerr << "[ParseEleIDVariables::initializeAnalyzer] " << DataEra << " is not implemented" << endl;
        exit(EXIT_FAILURE);
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, MCSample, IsDATA);
}

void ParseMuIDVariables::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> jets = GetAllJets();
    if (!PassNoiseFilter(jets, ev)) return;

    RVec<Electron> electrons = GetElectrons("POGTight", 25., 2.5);
    RVec<Muon> muons = GetMuons("", 10., 2.4);
    RVec<Gen> truth = GetAllGens();
    RVec<TrigObj> trigObjs = GetAllTrigObjs();

    if (! ev.PassTrigger(SglElTriggers)) return;
    if (! (electrons.size() == 1)) return;
    const auto &el = electrons.at(0);
    const float safePtCut = (Run == 3 ? 32: (DataYear == 2016 ? 30 : 35));
    if (! (el.Pt() > safePtCut)) return;
    if (! PassSLT(el, trigObjs)) return;
    if (! (muons.size() > 0)) return;
    
    // Update branches
    genWeight = MCweight()*ev.GetTriggerLumi("Full")*GetL1PrefireWeight()*myCorr->GetPUWeight(ev.nTrueInt());
    nMuons = muons.size();
    for (int i = 0; i < nMuons; i++) {
        const auto &mu = muons.at(i);
        pt[i] = mu.Pt();
        eta[i] = mu.Eta();
        isPOGMediumId[i] = mu.PassID("POGMedium");
        dZ[i] = mu.dZ();
        sip3d[i] = mu.SIP3D();
        tkRelIso[i] = mu.TkRelIso();
        miniPFRelIso[i] = mu.MiniPFRelIso();
        lepType[i] = GetLeptonType(mu, truth);

        // Use jetIdx for efficient jet matching
        nearestJetFlavour[i] = -1; // Default: no jet match
        
        short jetIdx = mu.JetIdx();
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

        isEMuTrigMatched[i] = PassEMT(mu, trigObjs);
        isIsoMuTrigMatched[i] = PassIsoMuT(mu, trigObjs);
    }
    Events->Fill();
}

void ParseMuIDVariables::WriteHist() {
    TFile* outfile = GetOutfile();
    Events->Write();
    outfile->Close();
}

bool ParseMuIDVariables::PassSLT(const Electron &el, const RVec<TrigObj> &trigObjs) {
    const float trig_pt_cut = (Run == 3 ? 30: (DataYear == 2016 ? 27 : 32));
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isElectron()) continue;
        if (! (trigObj.DeltaR(el) < 0.3)) continue;
        if (! (trigObj.hasBit(1))) continue;
        if (! (trigObj.Pt() > trig_pt_cut)) continue;
        return true;
    }
    return false;
}

bool ParseMuIDVariables::PassEMT(const Muon &mu, const RVec<TrigObj> &trigObjs) {
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(5))) continue;
        //if (! (trigObj.Pt() > pt_cut)) continue;
        return true;
    }
    return false;
}

bool ParseMuIDVariables::PassIsoMuT(const Muon &mu, const RVec<TrigObj> &trigObjs) {
    for (const auto &trigObj : trigObjs) {
        if (! trigObj.isMuon()) continue;
        if (! (trigObj.DeltaR(mu) < 0.3)) continue;
        if (! (trigObj.hasBit(0)))continue;
        return true;
    }
    return false;
}
