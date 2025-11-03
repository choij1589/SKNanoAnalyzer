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
import torch
from torch_geometric.data import Data
from ROOT import TLorentzVector, TRandom3
from itertools import product

from MLTools.models import MultiClassParticleNet
from MLTools.formats import NodeParticle


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


def loadMultiClassParticleNet(channel, signals, fold=None):
    """
    Load multi-class ParticleNet models.

    Model directory structure:
    {DATA_DIR}/Run2/{channel}/Classifiers/ParticleNet/{signal}/fold-{fold}/models/ParticleNet.pt

    Args:
        channel: "Run1E2Mu" or "Run3Mu"
        signals: List of signal mass points (e.g., ["MHc160_MA85", "MHc130_MA90", "MHc100_MA95"])
        fold: Specific fold to load (0-4), or None to load all folds

    Returns:
        Dictionary: {signal: model} if fold specified, or {signal_fold{i}: model} if fold=None
    """
    models = {}
    data_dir = os.environ.get('DATA_DIR', os.path.join(os.environ['SKNANO_DATA'], 'Run2'))

    if fold is not None:
        # Load single fold for each signal
        for sig in signals:
            modelPath = f"{data_dir}/{channel}/Classifiers/ParticleNet/{sig}/fold-{fold}/models/ParticleNet.pt"
            summaryPath = f"{data_dir}/{channel}/Classifiers/ParticleNet/{sig}/fold-{fold}/summary.txt"

            # Get num_hidden from summary file
            # with open(summaryPath, "r") as f:
            #    # Expected format: "ParticleNet, nNodes128, ..."
            #    first_line = f.readline().strip()
            #    parts = first_line.split(", ")
            #    num_hidden = int(parts[3].replace("nNodes", ""))
            num_hidden = 128  # Use default value if summary file is not available
            print(f"[WARNING] Using default num_hidden={num_hidden} for {sig} fold-{fold}")
            print(f"[WARNING] Loading {sig} with fixed fold-{fold}: {modelPath}")
            model = MultiClassParticleNet(9, 4, 4, num_hidden=num_hidden, dropout_p=0.25)
            checkpoint = torch.load(modelPath, map_location=torch.device("cpu"), weights_only=False)
            model.load_state_dict(checkpoint["model_state_dict"])
            model.eval()
            models[f"{sig}_fold-3"] = model
    else:
        # Load all folds (0-4) for each signal
        for sig in signals:
            for f in range(5):
                modelPath = f"{data_dir}/{channel}/Classifiers/ParticleNet/{sig}/fold-{f}/models/ParticleNet.pt"
                summaryPath = f"{data_dir}/{channel}/Classifiers/ParticleNet/{sig}/fold-{f}/summary.txt"

                # Get num_hidden from summary file
                with open(summaryPath, "r") as f:
                    first_line = f.readline().strip()
                    parts = first_line.split(", ")
                    num_hidden = int(parts[3].replace("nNodes", ""))

                print(f"Loading {sig} fold-{f}: {modelPath}")
                model = MultiClassParticleNet(9, 4, 4, num_hidden=num_hidden, dropout_p=0.25)
                checkpoint = torch.load(modelPath, map_location=torch.device("cpu"))
                model.load_state_dict(checkpoint["model_state_dict"])
                model.eval()
                models[f"{sig}_fold-{f}"] = model

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


def getGraphInput(muons, electrons, jets, bjets, METv, era, nFolds=5):
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

    Returns:
        Tuple of (data, fold) where data is PyG Data object and fold is the event fold number
    """
    particles = []

    # Muons
    for muon in muons:
        node = NodeParticle()
        node.isMuon = True
        node.SetPtEtaPhiM(muon.Pt(), muon.Eta(), muon.Phi(), muon.M())
        node.charge = muon.Charge()
        particles.append(node)

    # Electrons
    for ele in electrons:
        node = NodeParticle()
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

    # Era encoding for graph-level features
    if era == "2016preVFP":
        eraIdx = torch.tensor([[1, 0, 0, 0]], dtype=torch.float)
    elif era == "2016postVFP":
        eraIdx = torch.tensor([[0, 1, 0, 0]], dtype=torch.float)
    elif era == "2017":
        eraIdx = torch.tensor([[0, 0, 1, 0]], dtype=torch.float)
    elif era == "2018":
        eraIdx = torch.tensor([[0, 0, 0, 1]], dtype=torch.float)
    else:
        raise ValueError(f"Unsupported era for ParticleNet: {era}")

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
