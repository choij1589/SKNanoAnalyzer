#!/usr/bin/env python
"""
Data format classes for ParticleNet inference.

Provides NodeParticle wrapper class for ROOT physics objects with particle type flags.
"""

from ROOT import TLorentzVector
from ROOT import TMath


class NodeParticle(TLorentzVector):
    """
    TLorentzVector wrapper with particle type flags for ParticleNet node features.

    Extends TLorentzVector with additional properties needed for graph neural network input:
    - Particle type flags (muon, electron, jet, b-jet)
    - Electric charge
    - B-tagging score (for jets)
    """
    def __init__(self):
        TLorentzVector.__init__(self)
        self.charge = 0
        self.btagScore = 0.
        self.isMuon = False
        self.isElectron = False
        self.isJet = False
        self.isBjet = False

    def Charge(self):
        return self.charge

    def BtagScore(self):
        return self.btagScore

    def MT(self, part):
        """Calculate transverse mass with another particle."""
        dPhi = self.DeltaPhi(part)
        return TMath.Sqrt(2*self.Pt()*part.Pt()*(1-TMath.Cos(dPhi)))

    def IsMuon(self):
        return self.isMuon

    def IsElectron(self):
        return self.isElectron

    def IsJet(self):
        return self.isJet

    def IsBjet(self):
        return self.isBjet
