#!/usr/bin/env python
"""
DEPRECATED: This module is deprecated in favor of PromptAnalyzer.py

For SR + Tree only mode (equivalent to PromptTreeProducer):
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "NoHistMode"])

For SR + Tree + Theory systematics:
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "RunTheoryUnc", "NoHistMode"])

For SR + Both modes (default in PromptAnalyzer):
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst"])

Migration: Replace `from PromptTreeProducer import PromptTreeProducer` with
           `from PromptAnalyzer import PromptAnalyzer as PromptTreeProducer`

This file will be removed in a future version.

---
PromptTreeProducer: Python-based tree producer with ParticleNet multi-class inference.

Produces output trees with:
- Event kinematics (masses, transverse masses)
- ParticleNet scores for multiple mass points and background classes
- Fold information
- Event weights with systematic variations

This analyzer serves as a Python validation tool to verify ParticleNet C++ implementation.
"""
import warnings
warnings.warn(
    "PromptTreeProducer is deprecated. Use PromptAnalyzer instead. "
    "For SR + Tree only mode: Userflags = ['Run1E2Mu', 'RunSyst', 'NoHistMode']",
    DeprecationWarning,
    stacklevel=2
)
import os
from ROOT import TString
from ROOT.VecOps import RVec
from ROOT import TriLeptonBase
from ROOT import MyCorrection; myVar = MyCorrection.variation
from ROOT import SystematicHelper
from ROOT import JetTagging
from ROOT import Event, Lepton, Muon, Electron, Jet
from ROOT import Gen, GenJet
from ROOT import TTree
from array import array
from enum import IntEnum

from MLTools.helpers import loadMultiClassParticleNet, getGraphInput, getMultiClassScore


class CutStage(IntEnum):
    Initial = 0
    NoiseFilter = 1
    VetoMap = 2
    LeptonSelection = 3
    Trigger = 4
    KinematicCuts = 5
    JetRequirements = 6
    ConversionFilter = 7
    SamplePatching = 8
    Final = 9


class PromptTreeProducer(TriLeptonBase):
    """
    Python-based tree producer with ParticleNet multi-class scores.

    Implements the same event selection and object definitions as C++ PromptTreeProducer
    but uses Python-based ParticleNet inference for easier debugging and validation.
    """

    def __init__(self):
        super().__init__()
        self.systHelper = None
        self.models = {}  # ParticleNet models: {signal_fold{i}: model}
        self.trees = {}
        self.mass1 = {}
        self.mass2 = {}
        self.MT1 = {}
        self.MT2 = {}
        self.scores = {}  # Nested dict: [syst][signal][class]
        self.fold = {}
        self.weight = {}
        self.theorySystematics = []
        self.hasTheoryWeights = False

    def initializePyAnalyzer(self):
        """Initialize analyzer: channel flags, systematics, models, and trees."""
        self.initializeAnalyzer()

        # Channel flags
        if not (self.Run1E2Mu or self.Run3Mu):
            raise ValueError("Run1E2Mu or Run3Mu must be set")
        if self.Run1E2Mu and self.Run3Mu:
            raise ValueError("Run1E2Mu and Run3Mu cannot be set at the same time")

        # Determine channel
        self.channel = "Run1E2Mu" if self.Run1E2Mu else "Run3Mu"

        # ParticleNet configuration
        self.signals = ["MHc160_MA85", "MHc130_MA90", "MHc130_MA100", "MHc100_MA95", "MHc115_MA87", "MHc145_MA92", "MHc160_MA98"]
        self.classNames = ["signal", "nonprompt", "diboson", "ttZ"]

        # Load ParticleNet models (all folds for cross-validation)
        print(f"\n[PromptTreeProducer] Loading ParticleNet models for {self.channel}")
        self.models = loadMultiClassParticleNet(self.signals)
        print(f"[PromptTreeProducer] Loaded {len(self.models)} models")

        # Systematics
        if self.IsDATA:
            self.systHelper = SystematicHelper(
                f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml",
                self.DataStream,
                self.DataEra
            )
        elif not self.RunSyst:
            self.systHelper = SystematicHelper(
                f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml",
                self.MCSample,
                self.DataEra
            )
        else:
            self.systHelper = SystematicHelper(
                f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/TriLeptonSystematics.yaml",
                self.MCSample,
                self.DataEra
            )

        # Initialize theory systematics if enabled
        if self.RunTheoryUnc and not self.IsDATA:
            self.initializeTheorySystematics()

        # Prepare output trees
        self.__prepareTrees()

    def executeEvent(self):
        """Main event processing loop."""
        ev = self.GetEvent()

        # Initial cutflow entry
        initialWeight = 1.0 if self.IsDATA else self.MCweight() * ev.GetTriggerLumi("Full")

        # Event filters
        rawJets = self.GetAllJets()
        if not self.PassNoiseFilter(rawJets, ev):
            return

        rawMuons = self.GetAllMuons()
        if not (self.RunNoJetVeto or self.PassVetoMap(rawJets, rawMuons, "jetvetomap")):
            return

        rawElectrons = self.GetAllElectrons()
        truth = self.GetAllGens()
        genJets = self.GetAllGenJets()

        # Check for theory weights availability
        self.hasTheoryWeights = self.HasTheoryWeights()

        # Cache base object copies once per event
        self._cached_objects = {}

        def processEvent(syst, apply_weight_variation=False):
            """Process single systematic variation."""
            recoObjects = self.defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, syst)
            channelResult = self.selectEvent(ev, recoObjects, truth, syst)
            if channelResult is None:
                return None, None, None

            if apply_weight_variation:
                assert syst == "Central", "Only Central weight variation is allowed"
                weights = self.getWeights(ev, recoObjects, genJets, "Central")
                self.fillTree(channelResult, recoObjects, weights, "Central")

                # Get weight-only systematics
                weightOnlySysts = list(self.systHelper.getWeightOnlySystematics())
                if self.RunSyst and "1E2Mu" in channelResult:
                    weightOnlySysts.remove("DblMuTrigSF")
                if self.RunSyst and "3Mu" in channelResult:
                    weightOnlySysts.remove("ElectronIDSF")
                    weightOnlySysts.remove("EMuTrigSF")

                for systName in weightOnlySysts:
                    weights_up = self.getWeights(ev, recoObjects, genJets, f"{systName}_Up")
                    weights_down = self.getWeights(ev, recoObjects, genJets, f"{systName}_Down")
                    self.fillTree(channelResult, recoObjects, weights_up, f"{systName}_Up")
                    self.fillTree(channelResult, recoObjects, weights_down, f"{systName}_Down")

                # Return objects and weights for theory systematics processing
                return channelResult, recoObjects, weights
            else:
                weights = self.getWeights(ev, recoObjects, genJets, syst)
                self.fillTree(channelResult, recoObjects, weights, syst)
                return channelResult, recoObjects, weights

        # Process Central with weight variations
        channelResult, centralObjects, centralWeights = processEvent("Central", apply_weight_variation=True)

        # Process theory systematics (if enabled and weights available)
        if channelResult is not None and self.RunTheoryUnc and self.hasTheoryWeights:
            self.processTheorySystematics(channelResult, centralObjects, centralWeights)

        # Process systematics requiring evtLoopAgain
        for syst in self.systHelper:
            systName = syst.iter_name

            # Skip Central (already processed) and weight-only systematics
            if systName == "Central":
                continue
            if not self.systHelper.findSystematic(syst.syst_name).evtLoopAgain:
                continue

            # Process the systematic variation
            processEvent(systName, apply_weight_variation=False)

        # Clear cache at end of event
        self._cached_objects.clear()

    def defineObjects(self, ev: Event,
                            rawMuons: RVec[Muon],
                            rawElectrons: RVec[Electron],
                            rawJets: RVec[Jet],
                            genJets: RVec[GenJet],
                            syst):
        """
        Define physics objects with systematic variations.

        Returns dict with all selected objects and METv.
        """
        # Use cached object copies if available, otherwise create and cache them
        if syst == "Central":
            self._cached_objects["allMuons"] = RVec(Muon)(rawMuons)
            self._cached_objects["allElectrons"] = RVec(Electron)(rawElectrons)
            self._cached_objects["allJets"] = RVec(Jet)(rawJets)
            allMuons = self._cached_objects["allMuons"]
            allElectrons = self._cached_objects["allElectrons"]
            allJets = self._cached_objects["allJets"]
        else:
            # For systematic variations, start with cached copies and apply variations
            allMuons = RVec(Muon)(self._cached_objects.get("allMuons", rawMuons))
            allElectrons = RVec(Electron)(self._cached_objects.get("allElectrons", rawElectrons))
            allJets = RVec(Jet)(self._cached_objects.get("allJets", rawJets))

        # Apply scale variations
        if "ElectronEn" in syst:
            variation = syst.split("_")[1].lower()  # up or down
            allElectrons = self.ScaleElectrons(ev, allElectrons, variation)
        if "ElectronRes" in syst:
            variation = syst.split("_")[1].lower()  # up or down
            allElectrons = self.SmearElectrons(allElectrons, variation)
        if "MuonEn" in syst:
            variation = syst.split("_")[1].lower()  # up or down
            allMuons = self.ScaleMuons(allMuons, variation)
        if "JetEn" in syst:
            variation = syst.split("_")[1].lower()  # up or down
            allJets = self.ScaleJets(allJets, variation, "total")
        if "JetRes" in syst:
            variation = syst.split("_")[1].lower()  # up or down
            allJets = self.SmearJets(allJets, genJets, variation)

        MET_var = Event.MET_Syst.CENTRAL
        if "UnclusteredEn" in syst:
            if syst.split("_")[1].lower() == "up":
                MET_var = Event.MET_Syst.UE_UP
            else:
                MET_var = Event.MET_Syst.UE_DOWN

        METv = ev.GetMETVector(Event.MET_Type.PUPPI, MET_var)
        METv = self.ApplyTypeICorrection(METv, allJets, allElectrons, allMuons)

        # Sort objects in pt order
        allMuons = RVec(Muon)(sorted(allMuons, key=lambda x: x.Pt(), reverse=True))
        allElectrons = RVec(Electron)(sorted(allElectrons, key=lambda x: x.Pt(), reverse=True))
        allJets = RVec(Jet)(sorted(allJets, key=lambda x: x.Pt(), reverse=True))

        # Select objects
        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        tightMuons = self.SelectMuons(vetoMuons, self.MuonIDs.GetID("tight"), 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs.GetID("loose"), 10., 2.5)
        tightElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs.GetID("tight"), 15., 2.5)

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
                if not (self.RunNoJetVeto or self.PassVetoMap(j, allMuons, "jetvetomap")):
                    continue
            jets.emplace_back(j)

            # B-tagging
            if j.GetBTaggerResult(tagger) > wp:
                bjets.emplace_back(j)

        return {
            "vetoMuons": vetoMuons,
            "tightMuons": tightMuons,
            "vetoElectrons": vetoElectrons,
            "tightElectrons": tightElectrons,
            "jets_vetoLep": jets_vetoLep,
            "jets": jets,
            "bjets": bjets,
            "METv": METv
        }

    def selectEvent(self, ev: Event,
                          recoObjects: dict,
                          truth: RVec[Gen],
                          syst: str = "Central") -> str:
        """
        Apply event selection criteria.

        Returns channel name ("SR1E2Mu" or "SR3Mu") or None if event doesn't pass.
        """
        vetoMuons = recoObjects["vetoMuons"]
        tightMuons = recoObjects["tightMuons"]
        vetoElectrons = recoObjects["vetoElectrons"]
        tightElectrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        is1E2Mu = (tightElectrons.size() == 1 and vetoElectrons.size() == 1 and
                   tightMuons.size() == 2 and vetoMuons.size() == 2)
        is3Mu = (tightMuons.size() == 3 and vetoMuons.size() == 3 and
                 tightElectrons.size() == 0 and vetoElectrons.size() == 0)

        if self.Run1E2Mu and not is1E2Mu: return
        if self.Run3Mu and not is3Mu: return

        # Conversion sample filtering
        if self.MCSample.Contains("DYJets") or self.MCSample.Contains("TTG"):
            convMuons = RVec(Muon)()
            convElectrons = RVec(Electron)()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
            if not (convMuons.size() + convElectrons.size()) > 0: return

        # 1E2Mu baseline
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers): return

            mu1, mu2 = tightMuons.at(0), tightMuons.at(1)
            ele = tightElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return
            if not mu1.Charge() + mu2.Charge() == 0: return
            pair = mu1 + mu2
            if not pair.M() > 12.: return
            if not (jets.size() >= 2): return
            if not (bjets.size() >= 1): return
            return "SR1E2Mu"

        # 3Mu baseline
        if self.Run3Mu:
            if not ev.PassTrigger(self.DblMuTriggers): return

            mu1, mu2, mu3 = tuple(tightMuons)
            if not mu1.Pt() > 20.: return
            if not mu2.Pt() > 10.: return
            if not mu3.Pt() > 10.: return
            if not abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1: return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
            if not pair1.M() > 12.: return
            if not pair2.M() > 12.: return

            if not (jets.size() >= 2): return
            if not (bjets.size() >= 1): return 
            return "SR3Mu"
        return

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

    def getWeights(self, ev: Event,
                         recoObjects: dict,
                         genJets: RVec[GenJet],
                         syst="Central") -> dict:
        """Calculate all event weights and scale factors."""
        if self.IsDATA:
            return {
                "genWeight": 1.,
                "prefireWeight": 1.,
                "pileupWeight": 1.,
                "muonRecoSF": 1.,
                "muonIDSF": 1.,
                "eleRecoSF": 1.,
                "eleIDSF": 1.,
                "trigSF": 1.,
                "pileupIDSF": 1.,
                "btagSF": 1.,
                "WZNjetsSF": 1.
            }

        genWeight = self.MCweight() * ev.GetTriggerLumi("Full")

        var = myVar.nom
        if "Up" in syst:
            var = myVar.up
        if "Down" in syst:
            var = myVar.down

        if "L1Prefire" in syst:
            prefireWeight = self.GetL1PrefireWeight(var)
        else:
            prefireWeight = self.GetL1PrefireWeight(myVar.nom)

        if "PileupReweight" in syst:
            pileupWeight = self.myCorr.GetPUWeight(ev.nTrueInt(), var)
        else:
            pileupWeight = self.myCorr.GetPUWeight(ev.nTrueInt(), myVar.nom)

        muonRecoSF = self.myCorr.GetMuonRECOSF(recoObjects["tightMuons"])
        if "MuonIDSF" in syst:
            muonIDSF = self.myCorr.GetMuonIDSF("TopHNT", recoObjects["tightMuons"], var)
        else:
            muonIDSF = self.myCorr.GetMuonIDSF("TopHNT", recoObjects["tightMuons"])

        eleRecoSF = self.myCorr.GetElectronRECOSF(recoObjects["tightElectrons"])
        if "ElectronIDSF" in syst:
            eleIDSF = self.myCorr.GetElectronIDSF("TopHNT", recoObjects["tightElectrons"], var)
        else:
            eleIDSF = self.myCorr.GetElectronIDSF("TopHNT", recoObjects["tightElectrons"])

        trigSF = 1.
        if self.Run1E2Mu:
            if "EMuTrigSF" in syst:
                trigSF = self.myCorr.GetEMuTriggerSF(recoObjects["tightElectrons"], recoObjects["tightMuons"], var)
            else:
                trigSF = self.myCorr.GetEMuTriggerSF(recoObjects["tightElectrons"], recoObjects["tightMuons"])
        elif self.Run3Mu:
            if "DblMuTrigSF" in syst:
                trigSF = self.myCorr.GetDblMuTriggerSF(recoObjects["tightMuons"], var)
            else:
                trigSF = self.myCorr.GetDblMuTriggerSF(recoObjects["tightMuons"])

        pileupIDSF = 1.
        if self.Run == 2:
            jets = recoObjects["jets_vetoLep"]
            matched_idx = self.GenJetMatching(jets, genJets, self.fixedGridRhoFastjetAll, 0.4, 10.)
            if "PileupJetIDSF" in syst:
                pileupIDSF = self.myCorr.GetPileupJetIDSF(jets, matched_idx, "loose", var)
            else:
                pileupIDSF = self.myCorr.GetPileupJetIDSF(jets, matched_idx, "loose", myVar.nom)

        btagSF = 1.
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = JetTagging.JetFlavTaggerWP.Medium
        method = JetTagging.JetTaggingSFMethod.mujets
        source = "central"
        if "HFcorr" in syst:
            source = "hf_corr"
        if "HFuncorr" in syst:
            source = "hf_uncorr"
        if "LFcorr" in syst:
            source = "lf_corr"
        if "LFuncorr" in syst:
            source = "lf_uncorr"
        btagSF = self.myCorr.GetBTaggingReweightMethod1a(recoObjects["jets"], tagger, wp, method, var, source)

        WZNjetsSF = 1.
        if (self.Run == 3 and (self.MCSample.Contains("WZTo3LNu") or self.MCSample.Contains("ZZTo4L"))) and (not self.RunNoWZSF):
            njets = float(recoObjects["jets"].size())
            WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "Central")
            if "WZNjetsSF" in syst:
                if "Up" in syst:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "total_up")
                elif "Down" in syst:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "total_down")

        return {
            "genWeight": genWeight,
            "prefireWeight": prefireWeight,
            "pileupWeight": pileupWeight,
            "muonRecoSF": muonRecoSF,
            "muonIDSF": muonIDSF,
            "eleRecoSF": eleRecoSF,
            "eleIDSF": eleIDSF,
            "trigSF": trigSF,
            "pileupIDSF": pileupIDSF,
            "btagSF": btagSF,
            "WZNjetsSF": WZNjetsSF
        }

    def fillTree(self, channel: str,
                       recoObjects: dict,
                       weights: dict,
                       syst="Central"):
        """Fill output tree with event information and ParticleNet scores."""
        muons = recoObjects["tightMuons"]
        electrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        # Calculate total weight
        genWeight = weights["genWeight"]
        prefireWeight = weights["prefireWeight"]
        pileupWeight = weights["pileupWeight"]
        muonRecoSF = weights["muonRecoSF"]
        muonIDSF = weights["muonIDSF"]
        eleRecoSF = weights["eleRecoSF"]
        eleIDSF = weights["eleIDSF"]
        trigSF = weights["trigSF"]
        pileupIDSF = weights["pileupIDSF"]
        btagSF = weights["btagSF"]
        WZNjetsSF = weights["WZNjetsSF"]
        totWeight = genWeight * prefireWeight * pileupWeight
        totWeight *= muonRecoSF * muonIDSF * eleRecoSF * eleIDSF
        totWeight *= trigSF * pileupIDSF * btagSF * WZNjetsSF

        # Fill kinematic variables
        if "1E2Mu" in channel:
            mu1, mu2 = muons.at(0), muons.at(1)
            pair = mu1 + mu2
            self.mass1[syst][0] = pair.M()
            self.mass2[syst][0] = -999.
            ele = electrons.at(0)
            self.MT1[syst][0] = (ele + METv).Mt()
            self.MT2[syst][0] = -999.
        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1 = mu_ss1 + mu_os
            pair2 = mu_ss2 + mu_os
            self.mass1[syst][0] = pair1.M()
            self.mass2[syst][0] = pair2.M()
            self.MT1[syst][0] = (mu_ss1 + METv).Mt()
            self.MT2[syst][0] = (mu_ss2 + METv).Mt()

        # ParticleNet inference (only for Central systematic to save time)
        if syst == "Central":
            # Construct graph input
            data, fold = getGraphInput(muons, electrons, jets, bjets, METv, str(self.DataEra))
            self.fold[syst][0] = fold

            # Run inference for each signal mass point
            for signal in self.signals:
                if signal not in self.models.keys():
                    print(f"[WARNING] Model {signal} not found!")
                    continue

                model = self.models[signal]
                probs = getMultiClassScore(model, data)

                # Store scores: [P(signal), P(nonprompt), P(diboson), P(ttZ)]
                self.scores[syst][signal]["signal"][0] = probs[0]
                self.scores[syst][signal]["nonprompt"][0] = probs[1]
                self.scores[syst][signal]["diboson"][0] = probs[2]
                self.scores[syst][signal]["ttZ"][0] = probs[3]
        else:
            # For other systematics, copy scores from Central
            self.fold[syst][0] = self.fold["Central"][0]
            for signal in self.signals:
                for cls in self.classNames:
                    self.scores[syst][signal][cls][0] = self.scores["Central"][signal][cls][0]

        # Fill weight
        self.weight[syst][0] = totWeight

        # Fill tree
        self.trees[syst].Fill()

    def __prepareTrees(self):
        """Prepare output TTrees for all systematic variations."""
        # Get list of all systematics
        systList = ["Central"]
        if not self.IsDATA and self.RunSyst:
            # Add weight-only systematics with _Up and _Down variations
            for systName in self.systHelper.getWeightOnlySystematics():
                systList.append(f"{systName}_Up")
                systList.append(f"{systName}_Down")

            # Add systematics requiring evtLoopAgain
            for syst in self.systHelper:
                if syst.iter_name != "Central" and self.systHelper.findSystematic(syst.syst_name).evtLoopAgain:
                    systList.append(syst.iter_name)

        # Add theory systematics if enabled
        if self.RunTheoryUnc and not self.IsDATA:
            systList.extend(self.theorySystematics)

        for syst in systList:
            # Create TTree
            thisTree = TTree(f"Events_{syst}", "")

            # Kinematic branches
            self.mass1[syst] = array("d", [0.])
            thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.])
            thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            self.MT1[syst] = array("d", [0.])
            thisTree.Branch("MT1", self.MT1[syst], "MT1/D")
            self.MT2[syst] = array("d", [0.])
            thisTree.Branch("MT2", self.MT2[syst], "MT2/D")

            # ParticleNet score branches
            self.scores[syst] = {}
            for signal in self.signals:
                self.scores[syst][signal] = {}
                for cls in self.classNames:
                    self.scores[syst][signal][cls] = array("d", [0.])
                    branchName = f"score_{signal}_{cls}"
                    thisTree.Branch(branchName, self.scores[syst][signal][cls], f"{branchName}/D")

            # Fold branch
            self.fold[syst] = array("i", [0])
            thisTree.Branch("fold", self.fold[syst], "fold/I")

            # Weight branch
            self.weight[syst] = array("d", [0.])
            thisTree.Branch("weight", self.weight[syst], "weight/D")

            thisTree.SetDirectory(0)
            self.trees[syst] = thisTree

    def WriteHist(self):
        """Write all trees to output file."""
        self.GetOutfile().cd()
        for syst, tree in self.trees.items():
            tree.Write()

    # ===================================================================
    # Theory Uncertainty Methods
    # ===================================================================

    def initializeTheorySystematics(self):
        """Initialize list of theory systematic names."""
        self.theorySystematics = []

        # PDF variations (100): PDF_0 to PDF_99
        for i in range(100):
            self.theorySystematics.append(f"PDF_{i}")

        # Scale variations (7): indices 0,1,2,3,4,6,8 (skip 5 and 7)
        for i in [0, 1, 2, 3, 4, 6, 8]:
            self.theorySystematics.append(f"Scale_{i}")

        # PS variations (4)
        self.theorySystematics.extend(["PS_ISRUp", "PS_FSRUp", "PS_ISRDown", "PS_FSRDown"])

        # AlphaS variations (2)
        self.theorySystematics.extend(["AlphaS_Up", "AlphaS_Down"])

    def getTheoryWeight(self, systName: str) -> float:
        """Get theory weight for a given systematic variation."""
        # PDF variations: LHEPdfWeight[1-100] for PDF_0 to PDF_99
        if systName.startswith("PDF_"):
            idx = int(systName[4:])
            return self.LHEPdfWeight[idx + 1] if (idx + 1) < self.nLHEPdfWeight else 1.0

        # Scale variations: LHEScaleWeight[i]
        if systName.startswith("Scale_"):
            idx = int(systName[6:])
            return self.LHEScaleWeight[idx] if idx < self.nLHEScaleWeight else 1.0

        # PS variations: PSWeight[0-3]
        if systName == "PS_ISRUp":
            return self.PSWeight[0] if self.nPSWeight > 0 else 1.0
        if systName == "PS_FSRUp":
            return self.PSWeight[1] if self.nPSWeight > 1 else 1.0
        if systName == "PS_ISRDown":
            return self.PSWeight[2] if self.nPSWeight > 2 else 1.0
        if systName == "PS_FSRDown":
            return self.PSWeight[3] if self.nPSWeight > 3 else 1.0

        # AlphaS variations: LHEPdfWeight[101-102]
        if systName == "AlphaS_Up":
            return self.LHEPdfWeight[102] if self.nLHEPdfWeight > 102 else 1.0
        if systName == "AlphaS_Down":
            return self.LHEPdfWeight[101] if self.nLHEPdfWeight > 101 else 1.0

        return 1.0

    def processTheorySystematics(self, channel: str, recoObjects: dict, weights: dict):
        """Process theory systematic variations."""
        muons = recoObjects["tightMuons"]
        electrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        # Calculate base weight (without theory variations)
        baseWeight = (weights["genWeight"] * weights["prefireWeight"] * weights["pileupWeight"] *
                      weights["muonRecoSF"] * weights["muonIDSF"] * weights["eleRecoSF"] *
                      weights["eleIDSF"] * weights["trigSF"] * weights["pileupIDSF"] *
                      weights["btagSF"] * weights["WZNjetsSF"])

        for systName in self.theorySystematics:
            theoryWeight = self.getTheoryWeight(systName)
            finalWeight = baseWeight * theoryWeight

            # Fill kinematic variables (same as Central)
            if "1E2Mu" in channel:
                mu1, mu2 = muons.at(0), muons.at(1)
                pair = mu1 + mu2
                self.mass1[systName][0] = pair.M()
                self.mass2[systName][0] = -999.
                ele = electrons.at(0)
                self.MT1[systName][0] = (ele + METv).Mt()
                self.MT2[systName][0] = -999.
            elif "3Mu" in channel:
                mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
                pair1 = mu_ss1 + mu_os
                pair2 = mu_ss2 + mu_os
                self.mass1[systName][0] = pair1.M()
                self.mass2[systName][0] = pair2.M()
                self.MT1[systName][0] = (mu_ss1 + METv).Mt()
                self.MT2[systName][0] = (mu_ss2 + METv).Mt()

            # Copy fold and scores from Central (objects don't change for theory systematics)
            self.fold[systName][0] = self.fold["Central"][0]
            for signal in self.signals:
                for cls in self.classNames:
                    self.scores[systName][signal][cls][0] = self.scores["Central"][signal][cls][0]

            # Set weight
            self.weight[systName][0] = finalWeight

            # Fill tree
            self.trees[systName].Fill()
