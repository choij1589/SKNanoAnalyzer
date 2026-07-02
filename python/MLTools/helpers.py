#!/usr/bin/env python
"""
Helper functions for ParticleNet multi-class inference.

Provides utilities for:
- Loading trained models
- Constructing graph inputs with b-jet separation
- Computing k-NN graphs with DeltaR distance
- Running multi-class inference
- Fold calculation
"""

import os
import json
import re
import torch
from torch_geometric.data import Data
from ROOT import TLorentzVector, TRandom3

from MLTools.models import MultiClassParticleNet
from MLTools.formats import NodeParticle

SUPPORTED_PARTICLENET_MD_PREFIXES = ("MHc100", "MHc115", "MHc130", "MHc145", "MHc160")


def _particleNetMDRoot():
    return os.path.join(os.environ["SKNANO_DATA"], "All", "Combined", "Classifiers", "ParticleNetMD")


def _signalSortKey(signal):
    match = re.fullmatch(r"MHc(\d+)_MA(\d+)", signal)
    if match is None:
        return (10**9, 10**9, signal)
    return (int(match.group(1)), int(match.group(2)), signal)


def resolveParticleNetMDSignals(userflags):
    """
    Resolve ParticleNetMD mass-prefix userflags into concrete signal model names.

    Only prefix flags like "MHc160" are valid. Exact model flags such as
    "MHc160_MA85" are intentionally rejected to keep the job-level memory
    selector simple and explicit.
    """
    requested_prefixes = []
    for flag in userflags:
        flag = str(flag)
        if not flag.startswith("MHc"):
            continue
        if flag not in SUPPORTED_PARTICLENET_MD_PREFIXES:
            supported = ", ".join(SUPPORTED_PARTICLENET_MD_PREFIXES)
            raise ValueError(
                f"Unsupported ParticleNetMD userflag '{flag}'. "
                f"Use one of these mass prefixes: {supported}"
            )
        if flag not in requested_prefixes:
            requested_prefixes.append(flag)

    model_root = _particleNetMDRoot()
    if not os.path.isdir(model_root):
        raise FileNotFoundError(f"ParticleNetMD model directory does not exist: {model_root}")

    if not requested_prefixes:
        return []

    available_signals = [
        entry for entry in os.listdir(model_root)
        if os.path.isdir(os.path.join(model_root, entry))
    ]

    resolved = []
    for prefix in requested_prefixes:
        matches = sorted(
            [
                signal for signal in available_signals
                if re.fullmatch(rf"{re.escape(prefix)}_MA\d+", signal)
            ],
            key=_signalSortKey
        )
        if not matches:
            raise FileNotFoundError(
                f"No ParticleNetMD models found for userflag '{prefix}' under {model_root}"
            )

        for signal in matches:
            best_model_dir = os.path.join(model_root, signal, "best_model")
            model_path = os.path.join(best_model_dir, "model.pt")
            model_info_path = os.path.join(best_model_dir, "model_info.json")
            if not os.path.isfile(model_path) or not os.path.isfile(model_info_path):
                raise FileNotFoundError(
                    f"ParticleNetMD model '{signal}' is incomplete. "
                    f"Expected both {model_path} and {model_info_path}"
                )
            if signal not in resolved:
                resolved.append(signal)

    return sorted(resolved, key=_signalSortKey)


def getEdgeIndices(nodeList, k=4):
    """
    Compute k-nearest neighbors graph using DeltaR distance.

    CRITICAL: Uses DeltaR (angular distance in eta-phi space), NOT Euclidean distance!
    This must match the training data preprocessing in Preprocess.py.

    Args:
        nodeList: List of node feature vectors [[E, Px, Py, Pz, ...], ...]
        k: Number of nearest neighbors (default: 4)

    Returns:
        Tuple of (edge_index, edge_attribute) tensors
    """
    edgeIndex = []
    edgeAttribute = []

    for i, node in enumerate(nodeList):
        distances = {}
        for j, other in enumerate(nodeList):
            if node is other:  # Skip same node
                continue

            # Construct TLorentzVector for DeltaR calculation
            thisNode = TLorentzVector()
            otherNode = TLorentzVector()
            thisNode.SetPxPyPzE(node[1], node[2], node[3], node[0])
            otherNode.SetPxPyPzE(other[1], other[2], other[3], other[0])

            # Use DeltaR distance (matches training)
            distances[j] = thisNode.DeltaR(otherNode)

        # Sort by distance and take k nearest neighbors
        distances = dict(sorted(distances.items(), key=lambda item: item[1]))
        for idx in list(distances.keys())[:k]:
            edgeIndex.append([i, idx])
            edgeAttribute.append([distances[idx]])

    return torch.tensor(edgeIndex, dtype=torch.long), torch.tensor(edgeAttribute, dtype=torch.float)


def evtToGraph(nodeList, y=None, k=4):
    """
    Convert event node list to PyTorch Geometric Data object.

    Args:
        nodeList: List of node feature vectors
        y: Label (optional, for training)
        k: Number of nearest neighbors

    Returns:
        PyG Data object with x, edge_index, edge_attribute, and batch tensor
    """
    x = torch.tensor(nodeList, dtype=torch.float)
    edgeIndex, edgeAttribute = getEdgeIndices(nodeList, k=k)

    # Create batch tensor for single event (all nodes belong to graph 0)
    batch = torch.zeros(x.size(0), dtype=torch.long)

    graph = Data(
        x=x,
        y=y,
        edge_index=edgeIndex.t().contiguous(),
        edge_attribute=edgeAttribute,
        batch=batch  # Add batch tensor for proper GNN processing
    )
    return graph

def loadMultiClassParticleNetMD(signals):
    """
    Load mass-decorrelated multi-class ParticleNets (trained with DisCo loss).

    Model directory structure:
    {SKNANO_DATA}/All/Combined/Classifiers/ParticleNetMD/{signal}/best_model/model.pt

    Args:
        signals: List of signal mass points (e.g., ["MHc100_MA95", "MHc130_MA90", "MHc160_MA85"])
    Returns:
        Dictionary: {signal: model}
    """
    models = {}
    model_root = _particleNetMDRoot()

    for sig in signals:
        model_info_path = os.path.join(model_root, sig, "best_model", "model_info.json")
        with open(model_info_path, "r") as f:
            hyperparameters = json.load(f)["hyperparameters"]
        num_node_features = hyperparameters["num_node_features"]
        num_graph_features = hyperparameters["num_graph_features"]
        num_classes = hyperparameters["num_classes"]
        num_hidden = hyperparameters["num_hidden"]
        dropout_p = hyperparameters["dropout_p"]
        print(
            f"[INFO] Loaded ParticleNetMD config for {sig}: "
            f"nodes={num_node_features}, graph={num_graph_features}, "
            f"classes={num_classes}, hidden={num_hidden}, dropout={dropout_p}"
        )

        modelPath = os.path.join(model_root, sig, "best_model", "model.pt")
        print(f"Loading {sig}: {modelPath}")
        model = MultiClassParticleNet(
            num_node_features,
            num_graph_features,
            num_classes,
            num_hidden=num_hidden,
            dropout_p=dropout_p
        )
        checkpoint = torch.load(modelPath, map_location=torch.device("cpu"), weights_only=False)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.eval()
        model.sknano_num_graph_features = num_graph_features
        models[sig] = model

    return models
    
def calculateFold(METv, nJets, nFolds=5):
    """
    Calculate fold number for k-fold cross-validation.

    Args:
        METv: MET particle (Particle or TLorentzVector)
        nJets: Number of jets in the event
        nFolds: Number of folds (default: 5)

    Returns:
        Fold number (0 to nFolds-1)
    """
    randGen = TRandom3()
    seed = int(METv.Pt()) + 1
    randGen.SetSeed(seed)

    fold = -999
    for _ in range(nJets):
        fold = randGen.Integer(nFolds)

    return fold


def getGraphInput(muons, electrons, jets, bjets, METv, era, nFolds=5, num_graph_features=8):
    """
    Construct graph input with b-jet separated node features (Mode 2: separate_bjets=True).

    Node feature format: [E, Px, Py, Pz, charge, IsMuon, IsElectron, IsJet, IsBjet]

    Particle type encoding:
    - Muons:       [E, Px, Py, Pz, charge, 1, 0, 0, 0]
    - Electrons:   [E, Px, Py, Pz, charge, 0, 1, 0, 0]
    - Non-b-jets:  [E, Px, Py, Pz, 0     , 0, 0, 1, 0]
    - B-jets:      [E, Px, Py, Pz, 0,      0, 0, 1, 1]
    - MET:         [E, Px, Py, Pz, 0,      0, 0, 0, 0]

    IMPORTANT: This function handles the case where bjets are a SUBSET of jets.
    It will automatically filter out b-jets from the jets collection to avoid duplicates.

    Args:
        muons: ROOT vector of Muon objects
        electrons: ROOT vector of Electron objects
        jets: ROOT vector of Jet objects (ALL jets, including b-tagged)
        bjets: ROOT vector of Jet objects (b-tagged jets, subset of jets)
        METv: MET particle
        era: Data-taking era for graph-level features
        nFolds: Number of folds for fold calculation
        num_graph_features: Number of era-encoding features expected by the model

    Returns:
        Tuple of (data, fold) where data is PyG Data object and fold is the event fold number
    """
    particles = []

    swapLepton = (len(electrons) == 2 and len(muons) == 1)
    # Muons
    for muon in muons:
        node = NodeParticle()
        if swapLepton:
            node.isElectron = True
        else:
            node.isMuon = True 
        node.SetPtEtaPhiM(muon.Pt(), muon.Eta(), muon.Phi(), muon.M())
        node.charge = muon.Charge()
        particles.append(node)

    # Electrons
    for ele in electrons:
        node = NodeParticle()
        if swapLepton:
            node.isMuon = True
        else:
            node.isElectron = True
        node.SetPtEtaPhiM(ele.Pt(), ele.Eta(), ele.Phi(), ele.M())
        node.charge = ele.Charge()
        particles.append(node)

    # Create a set of b-jet identifiers (pt, eta, phi) to filter them out from jets
    bjet_ids = set()
    for bjet in bjets:
        # Use rounded values to handle floating point precision
        bjet_ids.add((round(bjet.Pt(), 6), round(bjet.Eta(), 6), round(bjet.Phi(), 6)))

    # Non-b-tagged jets (filter out bjets from jets collection)
    for jet in jets:
        jet_id = (round(jet.Pt(), 6), round(jet.Eta(), 6), round(jet.Phi(), 6))

        node = NodeParticle()
        node.isJet = True
        node.isBjet = jet_id in bjet_ids
        node.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M())
        node.charge = jet.Charge()
        particles.append(node)

    # MET as a node
    missing = NodeParticle()
    missing.SetPtEtaPhiM(METv.Pt(), 0., METv.Phi(), 0.)
    particles.append(missing)

    # Convert to node feature list
    nodeList = []
    for particle in particles:
        nodeList.append([
            particle.E(), particle.Px(), particle.Py(), particle.Pz(),
            particle.Charge(),
            1.0 if particle.IsMuon() else 0.0,
            1.0 if particle.IsElectron() else 0.0,
            1.0 if particle.IsJet() else 0.0,
            1.0 if particle.IsBjet() else 0.0
        ])
    
    # Create PyG Data object
    data = evtToGraph(nodeList, y=None, k=4)

    # Era encoding for graph-level features. ParticleNetMD uses all 8 eras.
    era_maps = {
        4: {"2016preVFP": 0, "2016postVFP": 1, "2017": 2, "2018": 3},
        8: {
            "2016preVFP": 0, "2016postVFP": 1, "2017": 2, "2018": 3,
            "2022": 4, "2022EE": 5, "2023": 6, "2023BPix": 7
        },
    }
    if num_graph_features not in era_maps:
        raise ValueError(f"Unsupported graph feature size {num_graph_features} for ParticleNet Input")
    era_map = era_maps[num_graph_features]
    if era not in era_map:
        raise ValueError(f"Unsupported era {era} for ParticleNet Input")
    eraIdx = torch.zeros(1, num_graph_features, dtype=torch.float)
    eraIdx[0, era_map[era]] = 1.0
    # Use consistent naming with training (graphInput, not graph_input)
    data.graphInput = eraIdx

    # Calculate fold
    # NOTE: Use len(jets) only, as bjets are a subset of jets in the inference code
    # This matches training where nJets represents total jet count
    fold = calculateFold(METv, len(jets), nFolds)

    return data, fold


def validateScoreDistribution(scores, tolerance=1e-5, verbose=False):
    """
    Validate that score distribution is proper probability distribution.

    Args:
        scores: NumPy array of class probabilities
        tolerance: Tolerance for sum check (default: 1e-5)
        verbose: Print validation details

    Returns:
        Boolean indicating if scores are valid

    Raises:
        ValueError if scores are invalid
    """
    import numpy as np

    # Check shape
    if scores.shape != (4,):
        raise ValueError(f"Score shape should be (4,), got {scores.shape}")

    # Check range [0, 1]
    if np.any(scores < 0) or np.any(scores > 1):
        raise ValueError(f"Scores out of [0,1] range! Min: {scores.min():.6f}, Max: {scores.max():.6f}")

    # Check sum to 1
    score_sum = np.sum(scores)
    if abs(score_sum - 1.0) > tolerance:
        raise ValueError(f"Scores don't sum to 1! Sum = {score_sum:.6f}")

    if verbose:
        print(f"[ParticleNet] Score validation passed:")
        print(f"  Sum: {score_sum:.6f}")
        print(f"  Range: [{scores.min():.6f}, {scores.max():.6f}]")
        print(f"  Distribution: Signal={scores[0]:.4f}, Nonprompt={scores[1]:.4f}, "
              f"Diboson={scores[2]:.4f}, ttZ={scores[3]:.4f}")

    return True


def getMultiClassScore(model, data, validate=True):
    """
    Run multi-class inference and return class probabilities.

    Args:
        model: Trained MultiClassParticleNet model
        data: PyG Data object with node features, edge_index, graphInput, and batch
        validate: Whether to validate score distribution (default: True)

    Returns:
        NumPy array of shape (4,) with class probabilities:
        [P(signal), P(nonprompt), P(diboson), P(ttZ)]
    """
    model.eval()
    with torch.no_grad():
        data.batch = torch.zeros(data.x.size(0), dtype=torch.long)

        # Forward pass with all required arguments (matching training)
        logits = model(data.x, data.edge_index, data.graphInput, data.batch)

        # Apply softmax to get probabilities
        probs = torch.softmax(logits, dim=1)
    scores = probs.numpy()[0]
    return scores
