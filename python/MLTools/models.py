#!/usr/bin/env python
"""
Multi-class ParticleNet models for 4-class classification.

Adapted from ChargedHiggsAnalysisV3/ParticleNet for SKNanoAnalyzer Python analysis.
Implements ParticleNet architecture for signal vs 3 background classification.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, LeakyReLU, Dropout, BatchNorm1d
from torch_geometric.nn import global_mean_pool, knn_graph
from torch_geometric.nn import GraphNorm
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import dropout_edge


class EdgeConv(MessagePassing):
    """
    EdgeConv layer for graph neural networks.

    Implements message passing between neighboring nodes with MLP transformation
    of concatenated node features.
    """
    def __init__(self, in_channels, out_channels, dropout_p):
        super().__init__(aggr="mean")
        self.mlp = Sequential(
            Linear(2*in_channels, out_channels),
            LeakyReLU(),
            BatchNorm1d(out_channels),
            Dropout(dropout_p),
            Linear(out_channels, out_channels),
            LeakyReLU(),
            BatchNorm1d(out_channels),
            Dropout(dropout_p),
            Linear(out_channels, out_channels),
            LeakyReLU(),
            BatchNorm1d(out_channels),
            Dropout(dropout_p)
        )

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, batch: torch.Tensor = None) -> torch.Tensor:
        return self.propagate(edge_index, x=x, batch=batch)

    def message(self, x_i: torch.Tensor, x_j: torch.Tensor) -> torch.Tensor:
        tmp = torch.cat([x_i, x_j - x_i], dim=1)
        return self.mlp(tmp)


class DynamicEdgeConv(EdgeConv):
    """
    Dynamic EdgeConv with k-nearest neighbors and residual connections.

    Builds k-NN graph dynamically and includes skip connections for better
    gradient flow. This is the core building block of the ParticleNet architecture.
    """
    def __init__(self, in_channels, out_channels, dropout_p, k=4):
        super().__init__(in_channels, out_channels, dropout_p=dropout_p)
        self.shortcut = Sequential(
            Linear(in_channels, out_channels),
            BatchNorm1d(out_channels),
            Dropout(dropout_p)
        )
        self.dropout_p = dropout_p
        self.k = k

    def forward(self, x, edge_index=None, batch=None):
        if edge_index is None:
            edge_index = knn_graph(x, self.k, batch, loop=False, flow=self.flow)
        edge_index, _ = dropout_edge(edge_index, p=self.dropout_p, training=self.training)
        out = super().forward(x, edge_index, batch=batch)
        out += self.shortcut(x)
        return out


class MultiClassParticleNet(torch.nn.Module):
    """
    Multi-class ParticleNet for 4-class classification.

    Architecture:
    - 3 DynamicEdgeConv layers with k=4 nearest neighbors
    - Concatenation of all conv outputs before pooling
    - Global mean pooling followed by dense layers
    - 4-class output (signal, nonprompt, diboson, ttZ)

    Input:
    - Node features: [E, Px, Py, Pz, charge, IsMuon, IsElectron, IsJet, IsBjet] (9 features)
    - Graph features: Era encoding (4 features)

    Output:
    - Logits for 4 classes (apply softmax to get probabilities)
    """
    def __init__(self, num_node_features, num_graph_features, num_classes=4, num_hidden=128, dropout_p=0.25):
        super(MultiClassParticleNet, self).__init__()

        # Input normalization
        self.gn0 = GraphNorm(num_node_features)

        # Three DynamicEdgeConv layers
        self.conv1 = DynamicEdgeConv(num_node_features, num_hidden, dropout_p, k=4)
        self.conv2 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, k=4)
        self.conv3 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, k=4)

        # Dense layers after global pooling
        self.bn0 = BatchNorm1d(num_hidden*3 + num_graph_features)
        self.dense1 = Linear(num_hidden*3 + num_graph_features, num_hidden)
        self.bn1 = BatchNorm1d(num_hidden)
        self.dense2 = Linear(num_hidden, num_hidden)
        self.bn2 = BatchNorm1d(num_hidden)
        self.output = Linear(num_hidden, num_classes)

        self.dropout_p = dropout_p
        self.num_classes = num_classes

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                graph_input: torch.Tensor, batch: torch.Tensor = None) -> torch.Tensor:
        """
        Forward pass of multi-class ParticleNet.

        Args:
            x: Node features (N_nodes, num_node_features)
            edge_index: Edge connectivity (2, N_edges)
            graph_input: Global graph features (N_graphs, num_graph_features)
            batch: Batch assignment for nodes (N_nodes,)

        Returns:
            Class logits (N_graphs, num_classes)
        """
        # Input normalization
        x = self.gn0(x, batch=batch)

        # Graph convolution layers with concatenation
        conv1 = self.conv1(x, edge_index, batch=batch)
        conv2 = self.conv2(conv1, batch=batch)
        conv3 = self.conv3(conv2, batch=batch)
        x = torch.cat([conv1, conv2, conv3], dim=1)

        # Global pooling
        x = global_mean_pool(x, batch=batch)
        x = torch.cat([x, graph_input], dim=1)
        x = self.bn0(x)

        # Dense classification layers
        x = F.leaky_relu(self.dense1(x))
        x = self.bn1(x)
        x = F.dropout(x, p=self.dropout_p, training=self.training)
        x = F.leaky_relu(self.dense2(x))
        x = self.bn2(x)
        x = F.dropout(x, p=self.dropout_p, training=self.training)
        x = self.output(x)

        # Return logits (softmax will be applied at inference time)
        return x
