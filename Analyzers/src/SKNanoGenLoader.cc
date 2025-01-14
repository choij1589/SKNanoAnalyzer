#include "SKNanoGenLoader.h"
using json = nlohmann::json;

SKNanoGenLoader::SKNanoGenLoader() {
    MaxEvent = -1;
    NSkipEvent = 0;
    LogEvery = 1000;
    MCSample = "";
    Campaign = "";
    xsec = 1.;
    sumW = 1.;
    sumSign = 1.;
    Userflags.clear();
}

SKNanoGenLoader::~SKNanoGenLoader() {
    delete nLHEPart;
    delete nGenPart;
    delete nGenVisTau;
    delete nGenJet;
    delete nGenJetAK8;
    delete nGenIsolatedPhoton;
    delete nGenDressedLepton;

    if (!fChain) return;
    if (fChain->GetCurrentFile()) fChain->GetCurrentFile()->Close();
}

void SKNanoGenLoader::Loop() {
    long nentries = fChain->GetEntries();
    if (MaxEvent > 0) nentries = std::min(nentries, MaxEvent);

    cout << "[SKNanoGenLoader::Loop] Event Loop Started" << endl;

    auto startTime = chrono::steady_clock::now();

    for (long jentry = 0; jentry < nentries; jentry++) {
        if (jentry < NSkipEvent) continue;

        if (jentry % LogEvery == 0) {
            auto currentTime = chrono::steady_clock::now();
            chrono::duration<double> elapsedTime = currentTime - startTime;
            double timePerEvent = elapsedTime.count() / (jentry + 1);
            double estimatedRemaining = (nentries - jentry) * timePerEvent;

            cout << "[SKNanoGenLoader::Loop] Processing " << jentry << " / " << nentries << " | Elapsed: " << fixed << setprecision(2) << elapsedTime.count() << "s" << ", Remaining: " << estimatedRemaining << "s" << endl;
        }
        if (fChain->GetEntry(jentry) < 0) exit(EIO);
        executeEvent();
    }

    cout << "[SKNanoGenLoader::Loop] Event Loop Finished" << endl;
}

void SKNanoGenLoader::SetMaxLeafSize() {
    df = new ROOT::RDataFrame(*fChain);

    auto getMaxBranchValue = [this](const TString &branchName) {
        if (!fChain->GetBranch(branchName)) {
            cout << "[SKNanoGenLoader::SetMaxLeafSize] Warning: Branch " << branchName << " not found" << endl;
            return 0;
        } 
        return static_cast<int>(*df->Max(branchName));
    };

    auto maxGenPartPtr = getMaxBranchValue("nGenPart");
    auto maxLHEPartPtr = getMaxBranchValue("nLHEPart");
    auto maxGenJetPtr = getMaxBranchValue("nGenJet");
    auto maxGenJetAK8Ptr = getMaxBranchValue("nGenJetAK8");
    auto maxGenIsolatedPhotonPtr = getMaxBranchValue("nGenIsolatedPhoton");
    auto maxGenDressedLeptonsPtr = getMaxBranchValue("nGenDressedLepton");
    auto maxGenVisTauPtr = getMaxBranchValue("nGenVisTau");
    
    // GenPart
    if (maxGenPartPtr > 0) {
        GenPart_eta.resize(maxGenPartPtr);
        GenPart_mass.resize(maxGenPartPtr);
        GenPart_pdgId.resize(maxGenPartPtr);
        GenPart_phi.resize(maxGenPartPtr);
        GenPart_pt.resize(maxGenPartPtr);
        GenPart_status.resize(maxGenPartPtr);
        GenPart_statusFlags.resize(maxGenPartPtr);
        GenPart_genPartIdxMother.resize(maxGenPartPtr);
    }

    // LHEPart
    if (maxLHEPartPtr > 0) {    
        LHEPart_pt.resize(maxLHEPartPtr);
        LHEPart_eta.resize(maxLHEPartPtr);
        LHEPart_phi.resize(maxLHEPartPtr);
        LHEPart_mass.resize(maxLHEPartPtr);
        LHEPart_incomingpz.resize(maxLHEPartPtr);
        LHEPart_pdgId.resize(maxLHEPartPtr);
        LHEPart_status.resize(maxLHEPartPtr);
        LHEPart_spin.resize(maxLHEPartPtr);
    }

    // GenJet
    if (maxGenJetPtr > 0) {
        GenJet_pt.resize(maxGenJetPtr);
        GenJet_eta.resize(maxGenJetPtr);
        GenJet_phi.resize(maxGenJetPtr);
        GenJet_mass.resize(maxGenJetPtr);
        GenJet_partonFlavour.resize(maxGenJetPtr);
        GenJet_hadronFlavour.resize(maxGenJetPtr);
    }

    // GenJetAK8
    if (maxGenJetAK8Ptr > 0) {
        GenJetAK8_pt.resize(maxGenJetAK8Ptr);
        GenJetAK8_eta.resize(maxGenJetAK8Ptr);
        GenJetAK8_phi.resize(maxGenJetAK8Ptr);
        GenJetAK8_mass.resize(maxGenJetAK8Ptr);
        GenJetAK8_partonFlavour.resize(maxGenJetAK8Ptr);
        GenJetAK8_hadronFlavour.resize(maxGenJetAK8Ptr);
    }

    // GenDressedLeptons
    if (maxGenDressedLeptonsPtr > 0) {
        GenDressedLepton_pt.resize(maxGenDressedLeptonsPtr);
        GenDressedLepton_eta.resize(maxGenDressedLeptonsPtr);
        GenDressedLepton_phi.resize(maxGenDressedLeptonsPtr);
        GenDressedLepton_mass.resize(maxGenDressedLeptonsPtr);
        GenDressedLepton_pdgId.resize(maxGenDressedLeptonsPtr);
        GenDressedLepton_hasTauAnc.resize(maxGenDressedLeptonsPtr);
    }

    // GenIsolatedPhoton
    if (maxGenIsolatedPhotonPtr > 0) {
        GenIsolatedPhoton_pt.resize(maxGenIsolatedPhotonPtr);
        GenIsolatedPhoton_eta.resize(maxGenIsolatedPhotonPtr);
        GenIsolatedPhoton_phi.resize(maxGenIsolatedPhotonPtr);
        GenIsolatedPhoton_mass.resize(maxGenIsolatedPhotonPtr);
    }

    // GenVisTau
    if (maxGenVisTauPtr > 0) {
        GenVisTau_pt.resize(maxGenVisTauPtr);
        GenVisTau_eta.resize(maxGenVisTauPtr);
        GenVisTau_phi.resize(maxGenVisTauPtr);
        GenVisTau_mass.resize(maxGenVisTauPtr);
        GenVisTau_charge.resize(maxGenVisTauPtr);
        GenVisTau_genPartIdxMother.resize(maxGenVisTauPtr);
        GenVisTau_status.resize(maxGenVisTauPtr);
    }

    delete df;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] Maximum Leaf Size Set" << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenPart: " << maxGenPartPtr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxLHEPart: " << maxLHEPartPtr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenJet: " << maxGenJetPtr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenJetAK8: " << maxGenJetAK8Ptr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenIsolatedPhoton: " << maxGenIsolatedPhotonPtr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenDressedLeptons: " << maxGenDressedLeptonsPtr << endl;
    cout << "[SKNanoGenLoader::SetMaxLeafSize] kMaxGenVisTau: " << maxGenVisTauPtr << endl;
}

void SKNanoGenLoader::Init() {
    cout << "[SKNanoGenLoader::Init] Initializing. Campaign = " << Campaign << endl;
    if(fChain->GetEntries() == 0){
        cout << "[SKNanoGenLoader::Init] No Entries in the Tree" << endl;
        cout << "[SKNanoGenLoader::Init] Exiting without make output..." << endl;
        exit(0);
    }
    SetMaxLeafSize();
    fChain->SetBranchStatus("*", 0);

    // Helper function to safely set branch address
    auto SafeSetBranchAddress = [this](const TString &branchName, void* address) {
        TBranch* branch = fChain->GetBranch(branchName);
        if (branch) {
            branch->SetAddress(address);
        } else {
            cout << "[SKNanoGenLoader::Init] Warning:Branch " << branchName << " not found" << endl;
        }
    };

    // Generator
    SafeSetBranchAddress("Generator_id1", &Generator_id1);
    SafeSetBranchAddress("Generator_id2", &Generator_id2);
    SafeSetBranchAddress("Generator_x1", &Generator_x1);
    SafeSetBranchAddress("Generator_x2", &Generator_x2);
    SafeSetBranchAddress("Generator_xpdf1", &Generator_xpdf1);
    SafeSetBranchAddress("Generator_xpdf2", &Generator_xpdf2);
    SafeSetBranchAddress("Generator_scalePDF", &Generator_scalePDF);
    SafeSetBranchAddress("Generator_weight", &Generator_weight);

    // LHE
    SafeSetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP);
    SafeSetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
    SafeSetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
    SafeSetBranchAddress("LHE_HT", &LHE_HT);
    SafeSetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming);
    SafeSetBranchAddress("LHE_Vpt", &LHE_Vpt);
    SafeSetBranchAddress("LHE_AlphaS", &LHE_AlphaS);
    SafeSetBranchAddress("LHE_Njets", &LHE_Njets);
    SafeSetBranchAddress("LHE_Nb", &LHE_Nb);
    SafeSetBranchAddress("LHE_Nc", &LHE_Nc);
    SafeSetBranchAddress("LHE_Nuds", &LHE_Nuds);
    SafeSetBranchAddress("LHE_Nglu", &LHE_Nglu);
    SafeSetBranchAddress("LHE_NpLO", &LHE_NpLO);
    SafeSetBranchAddress("LHE_NpNLO", &LHE_NpNLO);

    // LHEPart
    SafeSetBranchAddress("LHEPart_pt", LHEPart_pt.data());
    SafeSetBranchAddress("LHEPart_eta", LHEPart_eta.data());
    SafeSetBranchAddress("LHEPart_phi", LHEPart_phi.data());
    SafeSetBranchAddress("LHEPart_mass", LHEPart_mass.data());
    SafeSetBranchAddress("LHEPart_incomingpz", LHEPart_incomingpz.data());
    SafeSetBranchAddress("LHEPart_pdgId", LHEPart_pdgId.data());
    SafeSetBranchAddress("LHEPart_status", LHEPart_status.data());
    SafeSetBranchAddress("LHEPart_spin", LHEPart_spin.data());

    // GenPart
    SafeSetBranchAddress("GenPart_eta", GenPart_eta.data());
    SafeSetBranchAddress("GenPart_mass", GenPart_mass.data());
    SafeSetBranchAddress("GenPart_pdgId", GenPart_pdgId.data());
    SafeSetBranchAddress("GenPart_phi", GenPart_phi.data());
    SafeSetBranchAddress("GenPart_pt", GenPart_pt.data());
    SafeSetBranchAddress("GenPart_status", GenPart_status.data());
    SafeSetBranchAddress("GenPart_statusFlags", GenPart_statusFlags.data());
    SafeSetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother.data());

    // GenJet
    SafeSetBranchAddress("GenJet_pt", GenJet_pt.data());
    SafeSetBranchAddress("GenJet_eta", GenJet_eta.data());
    SafeSetBranchAddress("GenJet_phi", GenJet_phi.data());
    SafeSetBranchAddress("GenJet_mass", GenJet_mass.data());
    SafeSetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour.data());
    SafeSetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour.data());

    // GenJetAK8    
    SafeSetBranchAddress("GenJetAK8_pt", GenJetAK8_pt.data());
    SafeSetBranchAddress("GenJetAK8_eta", GenJetAK8_eta.data());
    SafeSetBranchAddress("GenJetAK8_phi", GenJetAK8_phi.data());
    SafeSetBranchAddress("GenJetAK8_mass", GenJetAK8_mass.data());
    SafeSetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour.data());
    SafeSetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour.data());

    // GenMET
    SafeSetBranchAddress("GenMET_pt", &GenMET_pt);
    SafeSetBranchAddress("GenMET_phi", &GenMET_phi);
    SafeSetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt);
    SafeSetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi);

    // GenDressedLeptons
    SafeSetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt.data());
    SafeSetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta.data());
    SafeSetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi.data());
    SafeSetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass.data());
    SafeSetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId.data());
    SafeSetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc.data());

    // GenIsolatedPhotons
    SafeSetBranchAddress("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt.data());
    SafeSetBranchAddress("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta.data());
    SafeSetBranchAddress("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi.data());
    SafeSetBranchAddress("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass.data());

    // GenVisTau
    SafeSetBranchAddress("GenVisTau_pt", GenVisTau_pt.data());
    SafeSetBranchAddress("GenVisTau_eta", GenVisTau_eta.data());
    SafeSetBranchAddress("GenVisTau_phi", GenVisTau_phi.data());
    SafeSetBranchAddress("GenVisTau_mass", GenVisTau_mass.data());
    SafeSetBranchAddress("GenVisTau_charge", GenVisTau_charge.data());
    SafeSetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother.data());
    SafeSetBranchAddress("GenVisTau_status", GenVisTau_status.data());
}