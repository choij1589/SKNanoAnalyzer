#include "GeneratorStudy.h"

GeneratorStudy::GeneratorStudy() {}
GeneratorStudy::~GeneratorStudy() {}

void GeneratorStudy::initializeAnalyzer() {
}
void GeneratorStudy::executeEvent() {
    LHEs = GetAllLHEs();
    GenJets = GetAllGenJets();
    GenDressedLeptons = GetAllGenDressedLeptons();

    // Split channels based on the number of outgoing jets in the LHE file
    RVec<LHE> lhe_leptons, lhe_muons, lhe_electrons, lhe_taus;
    RVec<LHE> lhe_outgoing_jets;
    for (const auto &lhe: LHEs) {
        if (lhe.Status() == -1) continue;

        if (fabs(lhe.PdgId()) == 11) {
            lhe_electrons.emplace_back(lhe);
            lhe_leptons.emplace_back(lhe);
        } else if (fabs(lhe.PdgId()) == 13) {
            lhe_muons.emplace_back(lhe);
            lhe_leptons.emplace_back(lhe);
        } else if (fabs(lhe.PdgId()) == 15) {
            lhe_taus.emplace_back(lhe);
            lhe_leptons.emplace_back(lhe);
        } else {
            lhe_outgoing_jets.emplace_back(lhe);
        }
    }
    
    TString n_outgoing_jets = "";
    if (lhe_outgoing_jets.size() == 0) {
        n_outgoing_jets = "0j";
    } else if (lhe_outgoing_jets.size() == 1) {
        n_outgoing_jets = "1j";
    } else if (lhe_outgoing_jets.size() == 2) {
        n_outgoing_jets = "2j";
    } else {
        cerr << "[GeneratorStudy::executeEvent] Unexpected number of outgoing jets" << endl;
    }

    TString channel = "DY"+n_outgoing_jets;
    if (lhe_electrons.size() == 2) {
        channel += "_ee";
    } else if (lhe_muons.size() == 2) {
        channel += "_mumu";
    } else if (lhe_taus.size() == 2) {
        channel += "_tautau";
    } else {
        cerr << "[GeneratorStudy::executeEvent] Unexpected lepton configuration" << endl;
    }
    const Particle lhe_ZCand = lhe_leptons[0] + lhe_leptons[1];

    // No matching for GenDressedLeptons
    RVec<GenDressedLepton> dressed_leptons, dressed_electrons, dressed_muons;
    if (GenDressedLeptons.size() == 2) {
        // Check at least one dressed lepton is from tau
        bool has_tau_anc = false;
        for (const auto &dressed_lep: GenDressedLeptons) {
            if (dressed_lep.HasTauAnc()) {
                has_tau_anc = true;
                break;
            }
        }

        if (!has_tau_anc) {
            for (const auto &dressed_lep: GenDressedLeptons) {
                dressed_leptons.emplace_back(dressed_lep);
                if (dressed_lep.PdgId() == 11) dressed_electrons.emplace_back(dressed_lep);
                else if (dressed_lep.PdgId() == 13) dressed_muons.emplace_back(dressed_lep);
            }
        }
    }
    const Particle dressed_ZCand = dressed_leptons[0] + dressed_leptons[1];
    
    // Fill LHE distributions
    const float weight = 1.0;
    FillHist("Inclusive/lhe/n_outgoing_jets", lhe_outgoing_jets.size(), weight, 5, 0., 5.);
    FillHist("Inclusive/lhe/Z_pt", lhe_ZCand.Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/Z_eta", lhe_ZCand.Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/Z_phi", lhe_ZCand.Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/Z_mass", lhe_ZCand.M(), weight, 100, 40., 140.);
    
    // LHE Lepton distributions
    FillHist("Inclusive/lhe/electrons/1/pt", lhe_electrons[0].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/electrons/1/eta", lhe_electrons[0].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/electrons/1/phi", lhe_electrons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/electrons/1/mass", lhe_electrons[0].M(), weight, 100, 40., 140.);
    FillHist("Inclusive/lhe/electrons/2/pt", lhe_electrons[1].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/electrons/2/eta", lhe_electrons[1].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/electrons/2/phi", lhe_electrons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/electrons/2/mass", lhe_electrons[1].M(), weight, 100, 40., 140.);
    
    FillHist("Inclusive/lhe/muons/1/pt", lhe_muons[0].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/muons/1/eta", lhe_muons[0].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/muons/1/phi", lhe_muons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/muons/1/mass", lhe_muons[0].M(), weight, 100, 40., 140.);
    FillHist("Inclusive/lhe/muons/2/pt", lhe_muons[1].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/muons/2/eta", lhe_muons[1].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/muons/2/phi", lhe_muons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/muons/2/mass", lhe_muons[1].M(), weight, 100, 40., 140.);

    FillHist("Inclusive/lhe/taus/1/pt", lhe_taus[0].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/taus/1/eta", lhe_taus[0].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/taus/1/phi", lhe_taus[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/taus/1/mass", lhe_taus[0].M(), weight, 100, 40., 140.);
    FillHist("Inclusive/lhe/taus/2/pt", lhe_taus[1].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/lhe/taus/2/eta", lhe_taus[1].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/lhe/taus/2/phi", lhe_taus[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/lhe/taus/2/mass", lhe_taus[1].M(), weight, 100, 40., 140.);

    // Fill dressed lepton distributions
    FillHist("Inclusive/gen/Z_pt", dressed_ZCand.Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/gen/Z_eta", dressed_ZCand.Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/gen/Z_phi", dressed_ZCand.Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/gen/Z_mass", dressed_ZCand.M(), weight, 100, 40., 140.);

    FillHist("Inclusive/gen/electrons/1/pt", dressed_electrons[0].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/gen/electrons/1/eta", dressed_electrons[0].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/gen/electrons/1/phi", dressed_electrons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/gen/electrons/1/mass", dressed_electrons[0].M(), weight, 100, 40., 140.);
    FillHist("Inclusive/gen/electrons/2/pt", dressed_electrons[1].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/gen/electrons/2/eta", dressed_electrons[1].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/gen/electrons/2/phi", dressed_electrons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/gen/electrons/2/mass", dressed_electrons[1].M(), weight, 100, 40., 140.);

    FillHist("Inclusive/gen/muons/1/pt", dressed_muons[0].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/gen/muons/1/eta", dressed_muons[0].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/gen/muons/1/phi", dressed_muons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/gen/muons/1/mass", dressed_muons[0].M(), weight, 100, 40., 140.);
    FillHist("Inclusive/gen/muons/2/pt", dressed_muons[1].Pt(), weight, 300, 0., 300.);
    FillHist("Inclusive/gen/muons/2/eta", dressed_muons[1].Eta(), weight, 100, -5., 5.);
    FillHist("Inclusive/gen/muons/2/phi", dressed_muons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist("Inclusive/gen/muons/2/mass", dressed_muons[1].M(), weight, 100, 40., 140.);

    FillHist(channel+"/lhe/Z_pt", lhe_ZCand.Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/Z_eta", lhe_ZCand.Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/Z_phi", lhe_ZCand.Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/Z_mass", lhe_ZCand.M(), weight, 100, 40., 140.);

    FillHist(channel+"/lhe/electrons/1/pt", lhe_electrons[0].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/electrons/1/eta", lhe_electrons[0].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/electrons/1/phi", lhe_electrons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/electrons/1/mass", lhe_electrons[0].M(), weight, 100, 40., 140.);
    FillHist(channel+"/lhe/electrons/2/pt", lhe_electrons[1].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/electrons/2/eta", lhe_electrons[1].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/electrons/2/phi", lhe_electrons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/electrons/2/mass", lhe_electrons[1].M(), weight, 100, 40., 140.);

    FillHist(channel+"/lhe/muons/1/pt", lhe_muons[0].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/muons/1/eta", lhe_muons[0].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/muons/1/phi", lhe_muons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/muons/1/mass", lhe_muons[0].M(), weight, 100, 40., 140.);
    FillHist(channel+"/lhe/muons/2/pt", lhe_muons[1].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/muons/2/eta", lhe_muons[1].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/muons/2/phi", lhe_muons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/muons/2/mass", lhe_muons[1].M(), weight, 100, 40., 140.);

    FillHist(channel+"/lhe/taus/1/pt", lhe_taus[0].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/taus/1/eta", lhe_taus[0].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/taus/1/phi", lhe_taus[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/taus/1/mass", lhe_taus[0].M(), weight, 100, 40., 140.);
    FillHist(channel+"/lhe/taus/2/pt", lhe_taus[1].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/lhe/taus/2/eta", lhe_taus[1].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/lhe/taus/2/phi", lhe_taus[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/lhe/taus/2/mass", lhe_taus[1].M(), weight, 100, 40., 140.);

    FillHist(channel+"/gen/Z_pt", dressed_ZCand.Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/gen/Z_eta", dressed_ZCand.Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/gen/Z_phi", dressed_ZCand.Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/gen/Z_mass", dressed_ZCand.M(), weight, 100, 40., 140.);
    
    FillHist(channel+"/gen/electrons/1/pt", dressed_electrons[0].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/gen/electrons/1/eta", dressed_electrons[0].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/gen/electrons/1/phi", dressed_electrons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/gen/electrons/1/mass", dressed_electrons[0].M(), weight, 100, 40., 140.);
    FillHist(channel+"/gen/electrons/2/pt", dressed_electrons[1].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/gen/electrons/2/eta", dressed_electrons[1].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/gen/electrons/2/phi", dressed_electrons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/gen/electrons/2/mass", dressed_electrons[1].M(), weight, 100, 40., 140.);

    FillHist(channel+"/gen/muons/1/pt", dressed_muons[0].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/gen/muons/1/eta", dressed_muons[0].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/gen/muons/1/phi", dressed_muons[0].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/gen/muons/1/mass", dressed_muons[0].M(), weight, 100, 40., 140.);
    FillHist(channel+"/gen/muons/2/pt", dressed_muons[1].Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/gen/muons/2/eta", dressed_muons[1].Eta(), weight, 100, -5., 5.);
    FillHist(channel+"/gen/muons/2/phi", dressed_muons[1].Phi(), weight, 64, -3.2, 3.2);
    FillHist(channel+"/gen/muons/2/mass", dressed_muons[1].M(), weight, 100, 40., 140.);
}