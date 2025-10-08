#include "TriLeptonBase.h"

TriLeptonBase::TriLeptonBase() : ModelPerFoldWarningIssued(false) {}
TriLeptonBase::~TriLeptonBase() {}

void TriLeptonBase::initializeAnalyzer() {
    // Flags
    Run1E2Mu = HasFlag("Run1E2Mu");
    Run3Mu = HasFlag("Run3Mu");
    RunNoVetoMap = HasFlag("RunNoVetoMap");
    RunNoWZSF = HasFlag("RunNoWZSF");
    RunSyst = HasFlag("RunSyst");

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

// ===================================================================
// ParticleNet GraphNet Implementation
// ===================================================================

void TriLeptonBase::loadGraphNetModels() {
    // Models were exported with torch.jit.trace() and expect edge_index as input
    // (computed in C++ via constructKNNGraph), so no torch_cluster dependency

    TString channel = Run1E2Mu ? "Run1E2Mu" : "Run3Mu";
    TString basePath = TString(getenv("SKNANO_DATA")) + "/Run2/" + channel + "/Classifiers/ParticleNet/";

    std::vector<TString> signals = {"MHc160_MA85", "MHc130_MA90", "MHc100_MA95"};

    for (const auto& sig : signals) {
        TString modelPath = basePath + sig + "_scripted.pt";
        try {
            GraphNetModels[sig] = torch::jit::load(modelPath.Data());
            GraphNetModels[sig].eval();
            std::cout << "[TriLeptonBase] Loaded ParticleNet model: " << sig << std::endl;
        } catch (const c10::Error& e) {
            std::cerr << "[ERROR] Failed to load ParticleNet model " << sig << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    // Issue warning about fold-specific models
    if (!ModelPerFoldWarningIssued) {
        std::cout << "\n[WARNING] ================================================" << std::endl;
        std::cout << "[WARNING] Currently loading single model per mass point." << std::endl;
        std::cout << "[WARNING] Future implementation should load separate models" << std::endl;
        std::cout << "[WARNING] for each fold (fold0-4) for proper k-fold inference." << std::endl;
        std::cout << "[WARNING] Expected path pattern: " << std::endl;
        std::cout << "[WARNING]   {signal}_fold{0-4}_scripted.pt" << std::endl;
        std::cout << "[WARNING] ================================================\n" << std::endl;
        ModelPerFoldWarningIssued = true;
    }
}

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

torch::Tensor TriLeptonBase::constructNodeFeatures(
    const RVec<Muon*>& muons,
    const RVec<Electron*>& electrons,
    const RVec<Jet*>& jets,
    const RVec<Jet*>& bjets,
    const Particle& METv)
{
    // Total particles: muons + electrons + jets + bjets + METv (1)
    int N = muons.size() + electrons.size() + jets.size() + bjets.size() + 1;

    auto x = torch::empty({N, 9}, torch::kFloat32);
    auto x_accessor = x.accessor<float, 2>();

    int idx = 0;

    // Muons: IsMuon=1, others=0
    for (const auto& mu : muons) {
        x_accessor[idx][0] = mu->E();
        x_accessor[idx][1] = mu->Px();
        x_accessor[idx][2] = mu->Py();
        x_accessor[idx][3] = mu->Pz();
        x_accessor[idx][4] = mu->Charge();
        x_accessor[idx][5] = 1.0;  // IsMuon
        x_accessor[idx][6] = 0.0;  // IsElectron
        x_accessor[idx][7] = 0.0;  // IsJet
        x_accessor[idx][8] = 0.0;  // IsBjet
        idx++;
    }

    // Electrons: IsElectron=1, others=0
    for (const auto& el : electrons) {
        x_accessor[idx][0] = el->E();
        x_accessor[idx][1] = el->Px();
        x_accessor[idx][2] = el->Py();
        x_accessor[idx][3] = el->Pz();
        x_accessor[idx][4] = el->Charge();
        x_accessor[idx][5] = 0.0;  // IsMuon
        x_accessor[idx][6] = 1.0;  // IsElectron
        x_accessor[idx][7] = 0.0;  // IsJet
        x_accessor[idx][8] = 0.0;  // IsBjet
        idx++;
    }

    // Jets (non-b-tagged): IsJet=1, IsBjet=0
    for (const auto& jet : jets) {
        x_accessor[idx][0] = jet->E();
        x_accessor[idx][1] = jet->Px();
        x_accessor[idx][2] = jet->Py();
        x_accessor[idx][3] = jet->Pz();
        x_accessor[idx][4] = jet->Charge();
        x_accessor[idx][5] = 0.0;  // IsMuon
        x_accessor[idx][6] = 0.0;  // IsElectron
        x_accessor[idx][7] = 1.0;  // IsJet = 1
        x_accessor[idx][8] = 0.0;  // IsBjet = 0
        idx++;
    }

    // BJets: IsJet=1, IsBjet=1 (BOTH TRUE)
    for (const auto& bjet : bjets) {
        x_accessor[idx][0] = bjet->E();
        x_accessor[idx][1] = bjet->Px();
        x_accessor[idx][2] = bjet->Py();
        x_accessor[idx][3] = bjet->Pz();
        x_accessor[idx][4] = bjet->Charge();
        x_accessor[idx][5] = 0.0;  // IsMuon
        x_accessor[idx][6] = 0.0;  // IsElectron
        x_accessor[idx][7] = 1.0;  // IsJet = 1 (TRUE for bjets)
        x_accessor[idx][8] = 1.0;  // IsBjet = 1 (TRUE for bjets)
        idx++;
    }

    // METv as a node: All flags = 0
    x_accessor[idx][0] = METv.E();
    x_accessor[idx][1] = METv.Px();
    x_accessor[idx][2] = METv.Py();
    x_accessor[idx][3] = METv.Pz();
    x_accessor[idx][4] = 0.0;  // Charge = 0
    x_accessor[idx][5] = 0.0;  // IsMuon = 0
    x_accessor[idx][6] = 0.0;  // IsElectron = 0
    x_accessor[idx][7] = 0.0;  // IsJet = 0
    x_accessor[idx][8] = 0.0;  // IsBjet = 0

    return x;
}

torch::Tensor TriLeptonBase::constructKNNGraph(const torch::Tensor& x, int k) {
    int N = x.size(0);
    std::vector<std::pair<int, int>> edges;

    auto x_accessor = x.accessor<float, 2>();

    // For each particle, find k nearest neighbors
    for (int i = 0; i < N; ++i) {
        std::vector<std::pair<float, int>> distances;

        for (int j = 0; j < N; ++j) {
            if (i != j) {
                // Euclidean distance in 9D feature space
                float dist_sq = 0.0;
                for (int d = 0; d < 9; ++d) {
                    float diff = x_accessor[i][d] - x_accessor[j][d];
                    dist_sq += diff * diff;
                }
                distances.push_back({dist_sq, j});
            }
        }

        // Sort and take k nearest
        std::sort(distances.begin(), distances.end());
        int num_neighbors = std::min(k, static_cast<int>(distances.size()));

        for (int n = 0; n < num_neighbors; ++n) {
            edges.push_back({i, distances[n].second});
        }
    }

    // Convert to edge_index tensor [2, N_edges]
    int num_edges = edges.size();
    auto edge_index = torch::empty({2, num_edges}, torch::kInt64);
    auto edge_accessor = edge_index.accessor<int64_t, 2>();

    for (int i = 0; i < num_edges; ++i) {
        edge_accessor[0][i] = edges[i].first;   // Source
        edge_accessor[1][i] = edges[i].second;  // Target
    }

    return edge_index;
}

torch::Tensor TriLeptonBase::getGraphFeatures(TString era) {
    // Era encoding: [2016preVFP, 2016postVFP, 2017, 2018]
    if (era == "2016preVFP") {
        return torch::tensor({{1.0, 0.0, 0.0, 0.0}});
    } else if (era == "2016postVFP") {
        return torch::tensor({{0.0, 1.0, 0.0, 0.0}});
    } else if (era == "2017") {
        return torch::tensor({{0.0, 0.0, 1.0, 0.0}});
    } else if (era == "2018") {
        return torch::tensor({{0.0, 0.0, 0.0, 1.0}});
    } else {
        std::cerr << "[ERROR] Unsupported era for ParticleNet: " << era << std::endl;
        std::cerr << "        Models trained on Run2 only (2016-2018)" << std::endl;
        exit(1);
    }
}

std::map<TString, std::vector<float>> TriLeptonBase::evalGraphNetScores(
    const RVec<Muon*>& muons,
    const RVec<Electron*>& electrons,
    const RVec<Jet*>& jets,
    const RVec<Jet*>& bjets,
    const Particle& METv,
    TString era)
{
    std::map<TString, std::vector<float>> scores;

    // Construct input tensors
    auto x = constructNodeFeatures(muons, electrons, jets, bjets, METv);
    auto edge_index = constructKNNGraph(x, 4);
    auto graph_input = getGraphFeatures(era);
    auto batch = torch::zeros({x.size(0)}, torch::kInt64);

    // Run inference for each signal model
    for (auto& [signal, model] : GraphNetModels) {
        std::vector<torch::jit::IValue> inputs = {x, edge_index, graph_input, batch};

        torch::Tensor output;
        {
            torch::NoGradGuard no_grad;  // Disable gradient computation
            output = model.forward(inputs).toTensor();
        }

        // Apply softmax to get probabilities
        auto probabilities = torch::softmax(output, 1);
        auto prob_accessor = probabilities.accessor<float, 2>();

        // Extract all 4 class probabilities: signal, nonprompt, diboson, ttZ
        std::vector<float> probs = {
            prob_accessor[0][0],  // signal
            prob_accessor[0][1],  // nonprompt
            prob_accessor[0][2],  // diboson
            prob_accessor[0][3]   // ttZ
        };
        scores[signal] = probs;
    }

    return scores;
}
