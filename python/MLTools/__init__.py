#!/usr/bin/env python
"""
MLTools package for ParticleNet multi-class inference.

Provides model definitions, data formats, and helper functions for
running ParticleNet inference in Python-based SKNanoAnalyzer analyzers.
"""

from MLTools.models import MultiClassParticleNet, EdgeConv, DynamicEdgeConv
from MLTools.formats import NodeParticle
from MLTools.helpers import (
    getEdgeIndices,
    evtToGraph,
    loadMultiClassParticleNet,
    calculateFold,
    getGraphInput,
    getMultiClassScore
)

__all__ = [
    # Models
    'MultiClassParticleNet',
    'EdgeConv',
    'DynamicEdgeConv',
    # Formats
    'NodeParticle',
    # Helpers
    'getEdgeIndices',
    'evtToGraph',
    'loadMultiClassParticleNet',
    'calculateFold',
    'getGraphInput',
    'getMultiClassScore',
]
