#include "Skim_TriLep.h"

Skim_TriLep::Skim_TriLep() {
    newtree = nullptr;
}

Skim_TriLep::~Skim_TriLep() {}

void Skim_TriLep::initializeAnalyzer() {
    GetOutfile()->cd();
    cout << "[Skim_TriLep::initializeAnalyzer] gDirectory = " << gDirectory->GetName() << endl;
    newtree = fChain->CloneTree(0);

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample, IsDATA);
}

void Skim_TriLep::executeEvent() {
    RVec<Muon> muonPreColl = GetMuons("NOCUT", 8., 2.4);
    RVec<Electron> electronPreColl = GetElectrons("NOCUT", 8., 2.4);

    // sort
    sort(muonPreColl.begin(), muonPreColl.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(electronPreColl.begin(), electronPreColl.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    
    int NMu = muonPreColl.size();
    int NEl = electronPreColl.size();
    int NLep = NEl+NMu;
    bool Has3Lep = (NLep >= 3);
    bool LeadLepPt15 = false;
    if (NMu > 0 && muonPreColl[0].Pt()>15) LeadLepPt15 = true;
    if (NEl > 0 && electronPreColl[0].Pt()>15) LeadLepPt15 = true;

    if (! (Has3Lep && LeadLepPt15) ) return;

    newtree->Fill();
}

void Skim_TriLep::WriteHist() {
    GetOutfile()->cd();
    newtree->Write();
}
