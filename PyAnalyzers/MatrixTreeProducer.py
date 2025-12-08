#!/usr/bin/env python
"""
MatrixTreeProducer: Python-based tree producer for matrix method fake background estimation.

Produces output tree with:
- Event kinematics (masses, transverse masses)
- ParticleNet scores for multiple mass points and background classes
- Fold information
- Fake weights from matrix method

Key differences from PromptTreeProducer:
- Uses LOOSE leptons (not tight) for matrix method
- Requires at least one NON-tight lepton per event
- NO systematic variations (Central only)
- Single output tree "Events"
- Uses GetFakeWeight instead of scale factors
"""
import os
from ROOT import TString
from ROOT.VecOps import RVec
from ROOT import TriLeptonBase
from ROOT import JetTagging
from ROOT import Event, Lepton, Muon, Electron, Jet
from ROOT import TTree
from array import array

from MLTools.helpers import loadMultiClassParticleNet, getGraphInput, getMultiClassScore


class MatrixTreeProducer(TriLeptonBase):
    """
    Python-based tree producer for matrix method with ParticleNet multi-class scores.

    Matrix method applies fake rates to events with loose (but not tight) leptons
    to estimate fake lepton backgrounds.
    """

    def __init__(self):
        super().__init__()
        self.models = {}  # ParticleNet models
        self.tree = None

        # Tree branch arrays
        self.run = array("i", [0])
        self.event = array("L", [0])
        self.lumi = array("i", [0])
        self.mass1 = array("d", [0.])
        self.mass2 = array("d", [0.])
        self.MT1 = array("d", [0.])
        self.MT2 = array("d", [0.])
        self.scores = {}  # Nested dict: [signal][class]
        self.fold = array("i", [0])
        self.weight = array("d", [0.])

    def initializePyAnalyzer(self):
        """Initialize analyzer: channel flags, models, and tree."""
        self.initializeAnalyzer()

        # Channel flags validation
        if not (self.Run1E2Mu or self.Run3Mu):
            raise ValueError("Run1E2Mu or Run3Mu must be set")
        if self.Run1E2Mu and self.Run3Mu:
            raise ValueError("Run1E2Mu and Run3Mu cannot be set at the same time")

        # Determine channel
        self.channel = "Run1E2Mu" if self.Run1E2Mu else "Run3Mu"

        # ParticleNet configuration
        self.signals = ["MHc160_MA85", "MHc130_MA90", "MHc130_MA100", "MHc100_MA95", "MHc115_MA87", "MHc145_MA92", "MHc160_MA98"]
        self.classNames = ["signal", "nonprompt", "diboson", "ttZ"]

        # Load ParticleNet models (fold=3 for testing)
        print(f"\n[MatrixTreeProducer] Loading ParticleNet models for {self.channel}")
        self.models = loadMultiClassParticleNet(self.signals)
        print(f"[MatrixTreeProducer] Loaded {len(self.models)} models")

        # Prepare output tree
        self.__prepareTree()

    def executeEvent(self):
        """Main event processing loop (no systematics)."""
        ev = self.GetEvent()

        # Event filters
        rawJets = self.GetAllJets()
        if not self.PassNoiseFilter(rawJets, ev):
            return

        rawMuons = self.GetAllMuons()
        if not (self.RunNoVetoMap or self.PassVetoMap(rawJets, rawMuons, "jetvetomap")):
            return

        rawElectrons = self.GetAllElectrons()

        # Define objects (no systematics, no genJets needed for matrix method)
        recoObjects = self.defineObjects(ev, rawMuons, rawElectrons, rawJets)

        # Select event
        channel = self.selectEvent(ev, recoObjects)

        if channel is not None:
            self.fillTree(channel, recoObjects)

    def defineObjects(self, ev: Event,
                            rawMuons: RVec[Muon],
                            rawElectrons: RVec[Electron],
                            rawJets: RVec[Jet]):
        """
        Define physics objects with veto, loose, and tight selections.

        Returns dict with veto/loose/tight leptons, jets, bjets, and METv.
        Matrix method requires all three selection levels!
        """
        # Create copies
        allMuons = RVec(Muon)(rawMuons)
        allElectrons = RVec(Electron)(rawElectrons)
        allJets = RVec(Jet)(rawJets)

        # MET correction
        METv = ev.GetMETVector(Event.MET_Type.PUPPI)
        METv = self.ApplyTypeICorrection(METv, allJets, allElectrons, allMuons)

        # Sort objects in pt order
        allMuons = RVec(Muon)(sorted(allMuons, key=lambda x: x.Pt(), reverse=True))
        allElectrons = RVec(Electron)(sorted(allElectrons, key=lambda x: x.Pt(), reverse=True))
        allJets = RVec(Jet)(sorted(allJets, key=lambda x: x.Pt(), reverse=True))

        # Select objects - CRITICAL: Keep veto, loose, AND tight selections
        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        looseMuons = self.SelectMuons(vetoMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        tightMuons = self.SelectMuons(looseMuons, self.MuonIDs.GetID("tight"), 10., 2.4)

        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs.GetID("loose"), 10., 2.5)
        looseElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs.GetID("loose"), 15., 2.5)
        tightElectrons = self.SelectElectrons(looseElectrons, self.ElectronIDs.GetID("tight"), 15., 2.5)

        max_jeteta = 2.4 if self.DataEra.Contains("2016") else 2.5
        jets_selected = self.SelectJets(allJets, "tight", 20., max_jeteta)
        jets_vetoLep = self.JetsVetoLeptonInside(jets_selected, vetoElectrons, vetoMuons, 0.4)

        # Single-pass filtering for jets and b-jets
        jets = RVec(Jet)()
        bjets = RVec(Jet)()
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = self.myCorr.GetBTaggingWP(tagger, JetTagging.JetFlavTaggerWP.Medium)

        for j in jets_vetoLep:
            # Apply Run 2 specific filters
            if self.Run == 2:
                if not j.PassID("loosePuId"):
                    continue
                if not (self.RunNoVetoMap or self.PassVetoMap(j, allMuons, "jetvetomap")):
                    continue
            jets.emplace_back(j)

            # B-tagging
            if j.GetBTaggerResult(tagger) > wp:
                bjets.emplace_back(j)

        return {
            "vetoMuons": vetoMuons,
            "looseMuons": looseMuons,
            "tightMuons": tightMuons,
            "vetoElectrons": vetoElectrons,
            "looseElectrons": looseElectrons,
            "tightElectrons": tightElectrons,
            "jets": jets,
            "bjets": bjets,
            "METv": METv
        }

    def selectEvent(self, ev: Event, recoObjects: dict) -> str:
        """
        Apply event selection criteria for matrix method.

        CRITICAL: Must have at least one NON-tight lepton for matrix method!
        Returns channel name ("SR1E2Mu" or "SR3Mu") or None.
        """
        vetoMuons = recoObjects["vetoMuons"]
        looseMuons = recoObjects["looseMuons"]
        tightMuons = recoObjects["tightMuons"]
        vetoElectrons = recoObjects["vetoElectrons"]
        looseElectrons = recoObjects["looseElectrons"]
        tightElectrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        # Check loose lepton multiplicity
        is1E2Mu = (looseElectrons.size() == 1 and vetoElectrons.size() == 1 and
                   looseMuons.size() == 2 and vetoMuons.size() == 2)
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)

        if self.Run1E2Mu and not is1E2Mu:
            return None
        if self.Run3Mu and not is3Mu:
            return None

        # CRITICAL: Require at least one non-tight lepton (matrix method requirement)
        if self.Run1E2Mu:
            if (tightMuons.size() == looseMuons.size() and
                tightElectrons.size() == looseElectrons.size()):
                return None

        if self.Run3Mu:
            if tightMuons.size() == looseMuons.size():
                return None

        # 1E2Mu baseline
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers):
                return None

            mu1, mu2 = looseMuons.at(0), looseMuons.at(1)
            ele = looseElectrons.at(0)

            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut:
                return None
            if not mu1.Charge() + mu2.Charge() == 0:
                return None

            pair = mu1 + mu2
            if not pair.M() > 12.:
                return None
            if not jets.size() >= 2:
                return None
            if not bjets.size() >= 1:
                return None

            return "SR1E2Mu"

        # 3Mu baseline
        if self.Run3Mu:
            if not ev.PassTrigger(self.DblMuTriggers):
                return None

            mu1, mu2, mu3 = tuple(looseMuons)
            if not mu1.Pt() > 20.:
                return None
            if not mu2.Pt() > 10.:
                return None
            if not mu3.Pt() > 10.:
                return None
            if not abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1:
                return None

            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(looseMuons)
            pair1 = mu_ss1 + mu_os
            pair2 = mu_ss2 + mu_os
            if not pair1.M() > 12.:
                return None
            if not pair2.M() > 12.:
                return None
            if not jets.size() >= 2:
                return None
            if not bjets.size() >= 1:
                return None

            return "SR3Mu"

        return None

    def configureChargeOf(self, muons: RVec[Muon]) -> tuple:
        """Configure muon charges for 3Mu channel."""
        if not muons.size() == 3:
            raise NotImplementedError(f"wrong no. of muons {muons.size()}")

        mu1, mu2, mu3 = tuple(muons)
        if mu1.Charge() == mu2.Charge():
            return (mu1, mu2, mu3)
        elif mu1.Charge() == mu3.Charge():
            return (mu1, mu3, mu2)
        elif mu2.Charge() == mu3.Charge():
            return (mu2, mu3, mu1)
        else:
            raise EOFError(f"wrong charge configuration {mu1.Charge()} {mu2.Charge()} {mu3.Charge()}")

    def fillTree(self, channel: str, recoObjects: dict):
        """Fill output tree with event information and ParticleNet scores."""
        looseMuons = recoObjects["looseMuons"]
        looseElectrons = recoObjects["looseElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        # Event identification
        self.run[0] = self.RunNumber
        self.event[0] = self.EventNumber
        self.lumi[0] = self.LumiBlock

        # Calculate fake weight (matrix method)
        self.weight[0] = self.GetFakeWeight(looseMuons, looseElectrons, "Central")

        # Fill kinematic variables using LOOSE leptons
        if "1E2Mu" in channel:
            mu1, mu2 = looseMuons.at(0), looseMuons.at(1)
            pair = mu1 + mu2
            self.mass1[0] = pair.M()
            self.mass2[0] = -999.
            self.MT1[0] = -999.
            self.MT2[0] = -999.
        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(looseMuons)
            pair1 = mu_ss1 + mu_os
            pair2 = mu_ss2 + mu_os
            self.mass1[0] = pair1.M()
            self.mass2[0] = pair2.M()
            MT1_val = (mu_ss1 + METv).Mt()
            MT2_val = (mu_ss2 + METv).Mt()
            self.MT1[0] = MT1_val
            self.MT2[0] = MT2_val

        # ParticleNet inference using LOOSE leptons
        data, fold = getGraphInput(looseMuons, looseElectrons, jets, bjets, METv, str(self.DataEra))
        self.fold[0] = fold

        # Run inference for each signal mass point
        for signal in self.signals:
            if signal not in self.models.keys():
                print(f"[WARNING] Model {signal} not found!")
                # Initialize scores to -999
                for cls in self.classNames:
                    self.scores[signal][cls][0] = -999.
                continue

            model = self.models[signal]
            probs = getMultiClassScore(model, data)

            # Store scores: [P(signal), P(nonprompt), P(diboson), P(ttZ)]
            self.scores[signal]["signal"][0] = probs[0]
            self.scores[signal]["nonprompt"][0] = probs[1]
            self.scores[signal]["diboson"][0] = probs[2]
            self.scores[signal]["ttZ"][0] = probs[3]

        # Fill tree
        self.tree.Fill()

    def __prepareTree(self):
        """Prepare single output TTree."""
        # Create TTree
        self.tree = TTree("Events", "Events")

        # Event identification
        self.tree.Branch("run", self.run, "run/I")
        self.tree.Branch("event", self.event, "event/l")
        self.tree.Branch("lumi", self.lumi, "lumi/I")

        # Kinematic branches
        self.tree.Branch("mass1", self.mass1, "mass1/D")
        self.tree.Branch("mass2", self.mass2, "mass2/D")
        self.tree.Branch("MT1", self.MT1, "MT1/D")
        self.tree.Branch("MT2", self.MT2, "MT2/D")

        # ParticleNet score branches
        for signal in self.signals:
            self.scores[signal] = {}
            for cls in self.classNames:
                self.scores[signal][cls] = array("d", [0.])
                branchName = f"score_{signal}_{cls}"
                self.tree.Branch(branchName, self.scores[signal][cls], f"{branchName}/D")

        # Fold and weight branches
        self.tree.Branch("fold", self.fold, "fold/I")
        self.tree.Branch("weight", self.weight, "weight/D")

        self.tree.SetDirectory(0)

    def WriteHist(self):
        self.GetOutfile().cd()
        self.tree.Write()