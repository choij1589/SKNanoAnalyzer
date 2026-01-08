#include "TriLeptonBase.h"

TriLeptonBase::TriLeptonBase() {}
TriLeptonBase::~TriLeptonBase() {}

void TriLeptonBase::initializeAnalyzer() {
    // Flags
    Run1E2Mu = HasFlag("Run1E2Mu");
    Run3Mu = HasFlag("Run3Mu");
    Run2E1Mu = HasFlag("Run2E1Mu");
    RunNoVetoMap = HasFlag("RunNoVetoMap");
    RunNoWZSF = HasFlag("RunNoWZSF");
    RunSyst = HasFlag("RunSyst");
    RunTheoryUnc = HasFlag("RunTheoryUnc");

    // Lepton IDs and triggers
    MuonIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    ElectronIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"
        };
    } else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
        };
    } else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", // prescaled
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            //"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8" // Need to measure the filter eff.
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
        };
    } else if (DataEra == "2018") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
        };
    } else {
       DblMuTriggers = {
           "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
       };
       EMuTriggers = {
           "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
           "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
       };
    }
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample, IsDATA);
}

void TriLeptonBase::executeEvent() {
    return;
}

float TriLeptonBase::GetFakeWeight(const RVec<Muon> &muons, const RVec<Electron> &electrons, const TString syst_key) {
    float weight = -1.;
    for (const auto &mu: muons) {
        const TString this_syst_key = syst_key.Contains("QCD") ? "QCD_MuEnriched" : syst_key;
        if (mu.PassID(MuonIDs->GetID("tight"))) continue;
        const float fr = myCorr->GetFakeRate(mu, "TopHNT", this_syst_key);
        weight *= -1.*(fr / (1.-fr));
    }

    for (const auto &ele: electrons) {
        const TString this_syst_key = syst_key.Contains("QCD") ? "QCD_EMEnriched" : syst_key;
        if (ele.PassID(ElectronIDs->GetID("tight"))) continue;
        const float fr = myCorr->GetFakeRate(ele, "TopHNT", this_syst_key);
        weight *= -1.*(fr / (1.-fr));
    }
    return weight;
}

RVec<Electron> TriLeptonBase::GetPTCorrScaledElectrons(const RVec<Electron> &electrons) {
    RVec<Electron> scaledElectrons;
    for (const auto &ele: electrons) {
        Electron scaledEle = ele;
        float ptCorr = ele.Pt()*(1.0+max(0., ele.MiniPFRelIso()-0.1));
        scaledEle.SetPtEtaPhiM(ptCorr, ele.Eta(), ele.Phi(), 0);
        scaledElectrons.emplace_back(scaledEle);
    }
    return scaledElectrons;
}

RVec<Muon> TriLeptonBase::GetPTCorrScaledMuons(const RVec<Muon> &muons) {
    RVec<Muon> scaledMuons;
    for (const auto &mu: muons) {
        Muon scaledMu = mu;
        float ptCorr = mu.Pt()*(1.0+max(0., mu.MiniPFRelIso()-0.1));
        scaledMu.SetPtEtaPhiM(ptCorr, mu.Eta(), mu.Phi(), mu.M());
        scaledMuons.emplace_back(scaledMu);
    }
    return scaledMuons;
}

// ===================================================================
// ParticleNet Helper Functions
// ===================================================================
int TriLeptonBase::calculateFold(const Particle& centralMETv, int nJets) {
    // Match Python implementation exactly:
    // randGen = TRandom3()
    // seed = int(METvPt) + 1
    // randGen.SetSeed(seed)
    // fold = -999
    // for _ in range(nJets):
    //     fold = randGen.Integer(nFolds)

    TRandom3 randGen;
    int seed = static_cast<int>(centralMETv.Pt()) + 1;  // +1 to avoid auto seed
    randGen.SetSeed(seed);

    int fold = -999;
    for (int i = 0; i < nJets; ++i) {
        fold = randGen.Integer(5);  // nFolds = 5
    }

    return fold;
}
