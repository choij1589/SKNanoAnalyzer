#!/usr/bin/env python
"""
PromptAnalyzer: Unified analyzer for prompt lepton analysis.

Supports multiple modes via Userflags:
- RunHistMode (default ON): Fill histograms
- RunTreeMode (default ON): Fill TTrees with ParticleNet scores
- NoHistMode: Disable histogram filling
- NoTreeMode: Disable tree filling
- RunCR: Use control region selection (WZ/ZG) instead of SR

Usage Examples:
    # SR + Both modes (default, recommended)
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst"])

    # SR + Histogram only
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "NoTreeMode"])

    # SR + Tree only
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "NoHistMode"])

    # CR + Both modes
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "RunCR"])

    # CR + Tree only
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "RunCR", "NoHistMode"])

    # Tree + theory systematics
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst", "RunTheoryUnc", "NoHistMode"])
"""
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
    EventVetoMap = 2
    LeptonSelection = 3
    ConversionFilter = 4
    Trigger = 5
    KinematicCuts = 6
    JetRequirements = 7
    JetVetoMap = 8
    Final = 9


class PromptAnalyzer(TriLeptonBase):
    """
    Unified analyzer for prompt lepton analysis supporting:
    - Signal region (SR) and control region (CR) selection
    - Histogram and TTree output modes
    - ParticleNet multi-class inference
    - Full systematic variations
    """

    def __init__(self):
        super().__init__()
        self.systHelper = None
        self.models = {}

        # Mode flags (set in _parseUserflags)
        self.RunHistMode = True   # Default: histogram mode
        self.RunTreeMode = False
        self.RunCR = False        # Default: SR mode (False = SR, True = CR)

        # Tree mode attributes
        self.trees = {}
        self.mass1 = {}
        self.mass2 = {}
        self.MT1 = {}
        self.MT2 = {}
        self.scores = {}
        self.fold = {}
        self.weight = {}
        self.theorySystematics = []
        self.hasTheoryWeights = False

    # ===================================================================
    # Initialization
    # ===================================================================

    def _parseUserflags(self):
        """Parse Userflags to set mode flags."""
        self.RunHistMode = True   # Default: both modes enabled
        self.RunTreeMode = True   # Default: both modes enabled
        self.RunCR = False        # Default: SR mode

        for flag in self.Userflags:
            flag_str = str(flag)
            if flag_str == "NoHistMode":
                self.RunHistMode = False
            if flag_str == "NoTreeMode":
                self.RunTreeMode = False
            if flag_str == "RunCR":
                self.RunCR = True

        if not self.RunHistMode and not self.RunTreeMode:
            raise ValueError("At least one output mode required (RunHistMode or RunTreeMode)")

        mode_str = []
        if self.RunHistMode:
            mode_str.append("Histogram")
        if self.RunTreeMode:
            mode_str.append("Tree")
        region_str = "CR (WZ/ZG)" if self.RunCR else "SR"
        print(f"[PromptAnalyzer] Mode: {'+'.join(mode_str)}, Region: {region_str}")

    def initializePyAnalyzer(self):
        """Initialize analyzer: channel flags, systematics, models, and trees."""
        self.initializeAnalyzer()

        # Parse mode flags
        self._parseUserflags()

        # Determine channel
        if self.Run1E2Mu:
            self.channel = "Run1E2Mu"
        elif self.Run3Mu:
            self.channel = "Run3Mu"
        elif self.Run2E1Mu:
            self.channel = "Run2E1Mu"
        else:
            raise ValueError("Run1E2Mu or Run3Mu or Run2E1Mu must be set")

        # ParticleNet configuration
        self.signals = ["MHc160_MA85", "MHc130_MA90", "MHc100_MA95", "MHc115_MA87", "MHc145_MA92", "MHc160_MA98"]
        self.classNames = ["signal", "nonprompt", "diboson", "ttZ"]

        # Load ParticleNet models
        print(f"[PromptAnalyzer] Loading ParticleNet models for {self.channel}")
        self.models = loadMultiClassParticleNet(self.signals)
        print(f"[PromptAnalyzer] Loaded {len(self.models)} models")

        # Systematics
        if self.IsDATA:
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml", self.DataStream, self.DataEra)
        elif not self.RunSyst:  # MC without RunSyst
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml", self.MCSample, self.DataEra)
        else:  # MC and RunSyst
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/TriLeptonSystematics.yaml", self.MCSample, self.DataEra)

        # Initialize theory systematics if enabled
        if self.RunTheoryUnc and not self.IsDATA:
            self.initializeTheorySystematics()

        # Prepare output trees if tree mode enabled
        if self.RunTreeMode:
            self.__prepareTrees()

    # ===================================================================
    # Cutflow and 2D Diagnostics (for histogram mode)
    # ===================================================================

    def fillCutflow(self, stage, channel, weight, syst):
        if not self.RunHistMode:
            return
        if not syst == "Central":
            return
        if weight is None:
            return
        cutIndex = int(stage)
        self.FillHist(f"{channel}/{syst}/cutflow", cutIndex, weight, 10, 0., 10.)

    def fillJetEtaPhi2D(self, jets, weight, stage):
        """Fill 2D jet eta-phi distribution for veto map validation."""
        if not self.RunHistMode:
            return
        for jet in jets:
            self.FillHist(f"JetEtaPhi/{stage}/eta_phi",
                          jet.Eta(), jet.Phi(), weight,
                          100, -5.0, 5.0,   # eta: 100 bins from -5 to 5
                          64, -3.2, 3.2)    # phi: 64 bins from -pi to pi

    def fillElectronScEtaPhi2D(self, electrons, weight, stage):
        """Fill 2D electron scEta-scPhi distribution for HEM veto validation."""
        if not self.RunHistMode:
            return
        for ele in electrons:
            self.FillHist(f"ElectronScEtaPhi/{stage}/scEta_scPhi",
                          ele.scEta(), ele.scPhi(), weight,
                          100, -2.5, 2.5,   # scEta: 100 bins from -2.5 to 2.5
                          64, -3.2, 3.2)    # scPhi: 64 bins from -pi to pi

    # ===================================================================
    # Main Event Loop
    # ===================================================================

    def executeEvent(self):
        ev = self.GetEvent()

        # Initial cutflow entry
        initialWeight = 1.0 if self.IsDATA else self.MCweight() * ev.GetTriggerLumi("Full")
        self.fillCutflow(CutStage.Initial, self.channel, initialWeight, "Central")

        rawJets = self.GetAllJets()
        if not self.PassNoiseFilter(rawJets, ev):
            return
        self.fillCutflow(CutStage.NoiseFilter, self.channel, initialWeight, "Central")

        rawMuons = self.GetAllMuons()
        if not (self.RunNoJetVeto or self.PassVetoMap(rawJets, rawMuons, "jetvetomap")):
            return
        self.fillCutflow(CutStage.EventVetoMap, self.channel, initialWeight, "Central")

        # Fill jet eta-phi for events passing veto map (Run 3)
        if self.Run == 3:
            self.fillJetEtaPhi2D(rawJets, initialWeight, "PassedEventVeto")

        rawElectrons = self.GetAllElectrons()
        truth = self.GetAllGens()
        genJets = self.GetAllGenJets()

        # Check for theory weights availability
        if self.RunTreeMode:
            self.hasTheoryWeights = self.HasTheoryWeights()

        # Cache base object copies once per event
        self._cached_objects = {}

        def processEvent(syst, apply_weight_variation=False):
            recoObjects = self.defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, syst)
            channel = self.selectEvent(ev, recoObjects, truth, syst, initialWeight if syst == "Central" else None)
            if channel is None:
                return None, None, None

            # Evaluate ParticleNet scores
            data, scores, fold = self.evalScore(
                recoObjects["tightMuons"],
                recoObjects["tightElectrons"],
                recoObjects["jets"],
                recoObjects["bjets"],
                recoObjects["METv"]
            )
            recoObjects["scores"] = scores
            recoObjects["fold"] = fold

            if apply_weight_variation:
                assert syst == "Central", "Only Central weight variation is allowed"
                weights = self.getWeights(ev, recoObjects, genJets, "Central")

                if self.RunHistMode:
                    self.fillObjects(channel, recoObjects, weights, "Central")
                if self.RunTreeMode:
                    self.fillTree(channel, recoObjects, weights, "Central")

                # Get weight-only systematics
                weightOnlySysts = list(self.systHelper.getWeightOnlySystematics())
                if self.RunSyst and "1E2Mu" in channel:
                    weightOnlySysts.remove("DblMuTrigSF")
                if self.RunSyst and "3Mu" in channel:
                    weightOnlySysts.remove("ElectronIDSF")
                    weightOnlySysts.remove("EMuTrigSF")

                for systName in weightOnlySysts:
                    weights_up = self.getWeights(ev, recoObjects, genJets, f"{systName}_Up")
                    weights_down = self.getWeights(ev, recoObjects, genJets, f"{systName}_Down")
                    if self.RunHistMode:
                        self.fillObjects(channel, recoObjects, weights_up, f"{systName}_Up")
                        self.fillObjects(channel, recoObjects, weights_down, f"{systName}_Down")
                    if self.RunTreeMode:
                        self.fillTree(channel, recoObjects, weights_up, f"{systName}_Up")
                        self.fillTree(channel, recoObjects, weights_down, f"{systName}_Down")

                return channel, recoObjects, weights
            else:
                weights = self.getWeights(ev, recoObjects, genJets, syst)
                if self.RunHistMode:
                    self.fillObjects(channel, recoObjects, weights, syst)
                if self.RunTreeMode:
                    self.fillTree(channel, recoObjects, weights, syst)
                return channel, recoObjects, weights

        # Process Central with weight variations
        channelResult, centralObjects, centralWeights = processEvent("Central", apply_weight_variation=True)

        # Process theory systematics (if enabled and tree mode)
        if channelResult is not None and self.RunTreeMode and self.RunTheoryUnc and self.hasTheoryWeights:
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

    # ===================================================================
    # Object Definition
    # ===================================================================

    def defineObjects(self, ev: Event,
                            rawMuons: RVec[Muon],
                            rawElectrons: RVec[Electron],
                            rawJets: RVec[Jet],
                            genJets: RVec[GenJet],
                            syst):
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
        # HEM veto: only apply to tight electrons, not veto (to properly reject events with extra leptons)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs.GetID("loose"), 10., 2.5)
        applyHEMVeto = (self.DataEra == "2018") and not self.RunNoHEMVeto

        # Fill electron scEta-scPhi BEFORE HEM veto (2018 only, Central only)
        if applyHEMVeto and syst == "Central":
            self.fillElectronScEtaPhi2D(vetoElectrons, 1.0, "BeforeHEMVeto")

        tightElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs.GetID("tight"), 15., 2.5, applyHEMVeto)

        # Fill electron scEta-scPhi AFTER HEM veto (2018 only, Central only)
        if applyHEMVeto and syst == "Central":
            self.fillElectronScEtaPhi2D(tightElectrons, 1.0, "AfterHEMVeto")

        max_jeteta = 2.4 if self.DataEra.Contains("2016") else 2.5
        jets_selected = self.SelectJets(allJets, "tight", 20., max_jeteta)
        jets_vetoLep = self.JetsVetoLeptonInside(jets_selected, vetoElectrons, vetoMuons, 0.4)
        jets_passPUID = RVec(Jet)()
        jets_vetoMap = RVec(Jet)()
        jets = RVec(Jet)()
        if self.Run == 2:
            jets_passPUID = self.SelectJets(jets_vetoLep, "loosePuId", 20., max_jeteta)

            # Fill jet eta-phi BEFORE jet-level veto (Run 2 only, Central only)
            if syst == "Central":
                self.fillJetEtaPhi2D(jets_passPUID, 1.0, "BeforeJetVeto")

            for j in jets_passPUID:
                if not (self.RunNoJetVeto or self.PassVetoMap(j, allMuons, "jetvetomap")):
                    continue
                jets_vetoMap.emplace_back(j)
            jets = jets_vetoMap

            # Fill jet eta-phi AFTER jet-level veto (Run 2 only, Central only)
            if syst == "Central":
                self.fillJetEtaPhi2D(jets, 1.0, "AfterJetVeto")
        else:   # Run3
            jets = jets_vetoLep

        # B-tagging
        bjets = RVec(Jet)()
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = self.myCorr.GetBTaggingWP(tagger, JetTagging.JetFlavTaggerWP.Medium)

        for j in jets:
            if j.GetBTaggerResult(tagger) > wp:
                bjets.emplace_back(j)

        return {"vetoMuons": vetoMuons,
                "tightMuons": tightMuons,
                "vetoElectrons": vetoElectrons,
                "tightElectrons": tightElectrons,
                "jets_vetoLep": jets_vetoLep,
                "jets_passPUID": jets_passPUID,
                "jets_vetoMap": jets_vetoMap,
                "jets": jets,
                "bjets": bjets,
                "METv": METv}

    # ===================================================================
    # Event Selection (Router)
    # ===================================================================

    def selectEvent(self, ev: Event,
                          recoObjects: dict,
                          truth: RVec[Gen],
                          syst: str = "Central",
                          weight: float = None) -> str:
        """Select event - routes to SR or CR selection based on RunCR flag."""
        if self.RunCR:
            return self._selectCREvent(ev, recoObjects, truth, syst, weight)
        else:
            return self._selectSREvent(ev, recoObjects, truth, syst, weight)

    def _selectSREvent(self, ev: Event,
                             recoObjects: dict,
                             truth: RVec[Gen],
                             syst: str = "Central",
                             weight: float = None) -> str:
        """Signal region event selection (from PromptSelector)."""
        vetoMuons = recoObjects["vetoMuons"]
        tightMuons = recoObjects["tightMuons"]
        vetoElectrons = recoObjects["vetoElectrons"]
        tightElectrons = recoObjects["tightElectrons"]
        jets_passPUID = recoObjects["jets_passPUID"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        is1E2Mu = (tightElectrons.size() == 1 and vetoElectrons.size() == 1 and
                   tightMuons.size() == 2 and vetoMuons.size() == 2)
        is3Mu = (tightMuons.size() == 3 and vetoMuons.size() == 3 and
                 tightElectrons.size() == 0 and vetoElectrons.size() == 0)
        is2E1Mu = (tightElectrons.size() == 2 and vetoElectrons.size() == 2 and
                   tightMuons.size() == 1 and vetoMuons.size() == 1)

        if self.Run1E2Mu and not is1E2Mu:
            return
        if self.Run3Mu and not is3Mu:
            return
        if self.Run2E1Mu and not is2E1Mu:
            return

        # Record lepton selection cutflow
        self.fillCutflow(CutStage.LeptonSelection, self.channel, weight, "Central")

        # For conversion samples
        if self.MCSample.Contains("DYJets") or self.MCSample.Contains("TTG") or self.MCSample.Contains("WWG"):
            # At least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = RVec(Muon)()
            convElectrons = RVec(Electron)()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) in [4, 5, -5, -6]:
                    convMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) in [4, 5, -5, -6]:
                    convElectrons.emplace_back(ele)
            if not (convMuons.size() + convElectrons.size()) > 0:
                return
        self.fillCutflow(CutStage.ConversionFilter, self.channel, weight, "Central")

        # 1E2Mu baseline
        # 1. pass EMuTriggers
        # 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        # 3. Exists OS muon pair with mass > 12 GeV
        # 4. At least two jets
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers):
                return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2 = tightMuons.at(0), tightMuons.at(1)
            ele = tightElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut:
                return
            if not mu1.Charge() + mu2.Charge() == 0:
                return
            pair = mu1 + mu2
            if not pair.M() > 12.:
                return

            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")
            if jets_passPUID.size() >= 2:
                self.fillCutflow(CutStage.JetRequirements, self.channel, weight, "Central")
            if not jets.size() >= 2:
                return
            self.fillCutflow(CutStage.JetVetoMap, self.channel, weight, "Central")

            if bjets.size() == 0:
                isOnZ = abs(pair.M() - 91.2) < 10.
                if isOnZ:
                    self.fillCutflow(CutStage.Final, "ZFake1E2Mu", weight, "Central")
                    return "ZFake1E2Mu"
                return
            else:
                self.fillCutflow(CutStage.Final, "SR1E2Mu", weight, "Central")
                return "SR1E2Mu"

        # 3Mu baseline
        # 1. pass DblMuTriggers
        # 2. Exact 3 tight muons, no additional leptons
        # 3. Exist OS muon pair,
        # 4. All OS muon pair mass > 12 GeV
        # 5. At least two jets
        if self.Run3Mu:
            if not ev.PassTrigger(self.DblMuTriggers):
                return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2, mu3 = tuple(tightMuons)
            if not mu1.Pt() > 20.:
                return
            if not mu2.Pt() > 10.:
                return
            if not mu3.Pt() > 10.:
                return
            if not abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1:
                return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
            if not pair1.M() > 12.:
                return
            if not pair2.M() > 12.:
                return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if jets_passPUID.size() >= 2:
                self.fillCutflow(CutStage.JetRequirements, self.channel, weight, "Central")
            if not jets.size() >= 2:
                return
            self.fillCutflow(CutStage.JetVetoMap, self.channel, weight, "Central")
            if bjets.size() == 0:
                isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
                if isOnZ:
                    self.fillCutflow(CutStage.Final, "ZFake3Mu", weight, "Central")
                    return "ZFake3Mu"
                return
            else:
                self.fillCutflow(CutStage.Final, "SR3Mu", weight, "Central")
                return "SR3Mu"

        # 2E1Mu sideband
        # 1. pass EMuTriggers
        # 2. Exact 2 tight electrons and 1 tight muon, no additional lepton
        # 3. Exist OS electron pair, 60 < M(ee) < 120 GeV
        # 4. At least two jets, at least one bjet
        if self.Run2E1Mu:
            if not ev.PassTrigger(self.EMuTriggers):
                return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            el1, el2 = tuple(tightElectrons)
            mu = tightMuons.at(0)
            passLeadMu = mu.Pt() > 25. and el1.Pt() > 15. and el2.Pt() > 15.
            passLeadEl = (el1.Pt() > 25. or el2.Pt() > 25.) and mu.Pt() > 10.
            passSafeCut = passLeadMu or passLeadEl
            if not passSafeCut:
                return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")
            if not (el1.Charge() + el2.Charge() == 0):
                return
            pair = el1 + el2
            if not 60. < pair.M() < 120.:
                return

            if not jets.size() >= 2:
                return
            if not bjets.size() >= 1:
                return
            self.fillCutflow(CutStage.Final, "TTZ2E1Mu", weight, "Central")
            return "TTZ2E1Mu"
        return

    def _selectCREvent(self, ev: Event,
                             recoObjects: dict,
                             truth: RVec[Gen],
                             syst: str = "Central",
                             weight: float = None) -> str:
        """Control region event selection (from CRPromptSelector)."""
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

        if self.Run1E2Mu and not is1E2Mu:
            return
        if self.Run3Mu and not is3Mu:
            return

        # Record lepton selection cutflow
        self.fillCutflow(CutStage.LeptonSelection, self.channel, weight, "Central")

        # For conversion samples
        if self.MCSample.Contains("DYJets") or self.MCSample.Contains("TTG"):
            # At least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = RVec(Muon)()
            convElectrons = RVec(Electron)()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) in [4, 5, -5, -6]:
                    convMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) in [4, 5, -5, -6]:
                    convElectrons.emplace_back(ele)
            if not (convMuons.size() + convElectrons.size()) > 0:
                return
            self.fillCutflow(CutStage.ConversionFilter, self.channel, weight, "Central")

        # 1E2Mu ZGamma
        # 1. pass EMuTriggers
        # 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        # 3. Exists OS muon pair with mass > 12 GeV
        # 4. |M(mumue) - 91.2| < 10
        # 5. No b-jet
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers):
                return
            if syst == "Central" and weight is not None:
                self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2 = tightMuons.at(0), tightMuons.at(1)
            ele = tightElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut:
                return
            if not mu1.Charge() + mu2.Charge() == 0:
                return
            pair = mu1 + mu2
            if not (pair.M() > 12.):
                return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair.M() - 91.2) < 10.):
                    return
                if not bjets.size() == 0:
                    return
                self.fillCutflow(CutStage.Final, "WZ1E2Mu", weight, "Central")
                return "WZ1E2Mu"
            else:
                if not (abs(pair.M() - 91.2) > 10.):
                    return
                if not abs((mu1 + mu2 + ele).M() - 91.2) < 10.:
                    return
                if not bjets.size() == 0:
                    return
                self.fillCutflow(CutStage.Final, "ZG1E2Mu", weight, "Central")
                return "ZG1E2Mu"

        # 3Mu ZGamma
        # 1. pass DblMuTriggers
        # 2. Exact 3 tight muons, no additional leptons
        # 3. Exist OS muon pair,
        # 4. All OS muon pair mass > 12 GeV
        # 5. |M(mumumu) - 91.2| < 10
        # 6. No b-jet
        if self.Run3Mu:
            if not ev.PassTrigger(self.DblMuTriggers):
                return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2, mu3 = tuple(tightMuons)
            if not mu1.Pt() > 20.:
                return
            if not mu2.Pt() > 10.:
                return
            if not mu3.Pt() > 10.:
                return
            if not abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1:
                return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
            if not pair1.M() > 12.:
                return
            if not pair2.M() > 12.:
                return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.):
                    return
                if not bjets.size() == 0:
                    return
                self.fillCutflow(CutStage.Final, "WZ3Mu", weight, "Central")
                return "WZ3Mu"
            else:
                if not abs(pair1.M() - 91.2) > 10.:
                    return
                if not abs(pair2.M() - 91.2) > 10.:
                    return
                if not abs((mu1 + mu2 + mu3).M() - 91.2) < 10.:
                    return
                if not bjets.size() == 0:
                    return
                self.fillCutflow(CutStage.Final, "ZG3Mu", weight, "Central")
                return "ZG3Mu"
        return

    # ===================================================================
    # Helper Methods
    # ===================================================================

    def configureChargeOf(self, muons: RVec[Muon]) -> tuple[Muon, Muon, Muon]:
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

    def getConvLepton(self, muons: RVec[Muon], electrons: RVec[Electron]) -> Lepton:
        """Get conversion lepton candidate for ZG CR (from CRPromptSelector)."""
        if electrons.size() == 1 and muons.size() == 2:
            return electrons.at(0)
        elif muons.size() == 3:
            mu1, mu2, mu3 = tuple(muons)
            if mu1.Charge() == mu2.Charge():
                pair1, pair2 = mu1 + mu3, mu2 + mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu2
            elif mu1.Charge() == mu3.Charge():
                pair1, pair2 = mu1 + mu2, mu2 + mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1, pair2 = mu1 + mu2, mu1 + mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu2
                else:
                    return mu3
        else:
            raise NotImplementedError(f"Wrong number of muons ({muons.size()})")

    def evalScore(self, muons, electrons, jets, bjets, METv):
        """Evaluate ParticleNet scores for all signal mass points."""
        scores = {}
        data, fold = getGraphInput(muons, electrons, jets, bjets, METv, str(self.DataEra))

        for sig in self.signals:
            if sig not in self.models.keys():
                print(f"[WARNING] Model {sig} not found!")
                for cls in self.classNames:
                    scores[f"{sig}_{cls}"] = -999.
                continue

            model = self.models[sig]
            probs = getMultiClassScore(model, data)

            # Store scores: [P(signal), P(nonprompt), P(diboson), P(ttZ)]
            scores[f"{sig}_signal"] = probs[0]
            scores[f"{sig}_nonprompt"] = probs[1]
            scores[f"{sig}_diboson"] = probs[2]
            scores[f"{sig}_ttZ"] = probs[3]

        return data, scores, fold

    # ===================================================================
    # Weight Calculation
    # ===================================================================

    def getWeights(self, ev: Event,
                         recoObjects: dict,
                         genJets: RVec[GenJet],
                         syst="Central") -> dict:
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
        if self.Run1E2Mu or self.Run2E1Mu:
            if "EMuTrigSF" in syst:
                trigSF = self.myCorr.GetEMuTriggerSF(recoObjects["tightElectrons"], recoObjects["tightMuons"], var)
            else:
                trigSF = self.myCorr.GetEMuTriggerSF(recoObjects["tightElectrons"], recoObjects["tightMuons"])
        elif self.Run3Mu:
            if "DblMuTrigSF" in syst:
                trigSF = self.myCorr.GetDblMuTriggerSF(recoObjects["tightMuons"], var)
            else:
                trigSF = self.myCorr.GetDblMuTriggerSF(recoObjects["tightMuons"])
        else:
            pass

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

        # WZ/ZZ NJets SF: apply to both WZTo3LNu and ZZTo4L
        WZNjetsSF = 1.
        if self.Run == 3 and self.MCSample.Contains("WZTo3LNu") and (not self.RunNoWZSF):
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

    # ===================================================================
    # Histogram Filling
    # ===================================================================

    def fillObjects(self, channel: str,
                          recoObjects: dict,
                          weights: dict,
                          syst="Central"):
        """Fill histograms - routes to SR or CR specific filling."""
        if not self.RunHistMode:
            return

        # Common histogram filling
        self._fillCommonObjects(channel, recoObjects, weights, syst)

        # Mode-specific histogram filling
        if self.RunCR:
            self._fillCRObjects(channel, recoObjects, weights, syst)
        else:
            self._fillSRObjects(channel, recoObjects, weights, syst)

    def _fillCommonObjects(self, channel: str,
                                 recoObjects: dict,
                                 weights: dict,
                                 syst="Central"):
        """Fill common histograms for both SR and CR."""
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
        totWeight *= trigSF * pileupIDSF * btagSF
        totWeight *= WZNjetsSF

        # Fill weights
        self.FillHist(f"{channel}/{syst}/weights/genWeight", genWeight, 1., 200, -10000, 10000.)
        self.FillHist(f"{channel}/{syst}/weights/prefireWeight", prefireWeight, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/pileupWeight", pileupWeight, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/muonRecoSF", muonRecoSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/muonIDSF", muonIDSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/eleRecoSF", eleRecoSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/eleIDSF", eleIDSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/trigSF", trigSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/pileupIDSF", pileupIDSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/btagSF", btagSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/WZNjetsSF", WZNjetsSF, 1., 100, -5., 5.)
        self.FillHist(f"{channel}/{syst}/weights/totWeight", totWeight, 1., 100, -5., 5.)

        # Fill base observables - muons
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), totWeight, 10, 0., 1.)

        # Fill base observables - electrons
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), totWeight, 100, 0., 1.)

        # Fill base observables - jets
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), totWeight, 100, 0., 100.)

        # Fill base observables - bjets
        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), totWeight, 100, 0., 100.)

        # Fill object counts
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), totWeight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), totWeight, 15, 0., 15.)

        # Fill MET
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), totWeight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), totWeight, 64, -3.2, 3.2)

    def _fillSRObjects(self, channel: str,
                             recoObjects: dict,
                             weights: dict,
                             syst="Central"):
        """Fill SR-specific histograms (from PromptSelector)."""
        muons = recoObjects["tightMuons"]
        electrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        # Calculate total weight
        totWeight = weights["genWeight"] * weights["prefireWeight"] * weights["pileupWeight"]
        totWeight *= weights["muonRecoSF"] * weights["muonIDSF"] * weights["eleRecoSF"] * weights["eleIDSF"]
        totWeight *= weights["trigSF"] * weights["pileupIDSF"] * weights["btagSF"]
        totWeight *= weights["WZNjetsSF"]

        # Fill additional SR-specific observables
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/energy", mu.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/px", mu.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/py", mu.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pz", mu.Pz(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/charge", mu.Charge(), totWeight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/muons/{idx}/min_dR_bjets", min([mu.DeltaR(bjet) for bjet in bjets]), totWeight, 100, 0., 10.)

        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/energy", ele.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/px", ele.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/py", ele.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pz", ele.Pz(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/charge", ele.Charge(), totWeight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/electrons/{idx}/min_dR_bjets", min([ele.DeltaR(bjet) for bjet in bjets]), totWeight, 100, 0., 10.)

        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/energy", jet.E(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/px", jet.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/py", jet.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pz", jet.Pz(), totWeight, 600, -300., 300.)

        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/energy", bjet.E(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/px", bjet.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/py", bjet.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pz", bjet.Pz(), totWeight, 600, -300., 300.)

        self.FillHist(f"{channel}/{syst}/METv/px", METv.Px(), totWeight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/py", METv.Py(), totWeight, 500, -250., 250.)

        # Fill pair masses and delta R
        if "1E2Mu" in channel:
            pair = muons.at(0) + muons.at(1)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), totWeight, 200, 0., 200.)
            # Delta R between leptons
            dR_ele_mu1 = electrons.at(0).DeltaR(muons.at(0))
            dR_ele_mu2 = electrons.at(0).DeltaR(muons.at(1))
            dR_mu1_mu2 = muons.at(0).DeltaR(muons.at(1))
            self.FillHist(f"{channel}/{syst}/dR_ele_mu1", dR_ele_mu1, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele_mu2", dR_ele_mu2, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_mu1_mu2", dR_mu1_mu2, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_min_ele_mu", min([dR_ele_mu1, dR_ele_mu2]), totWeight, 100, 0., 10.)

            # Within Z mass window?
            if 60 < pair.M() and pair.M() < 120:
                self.FillHist(f"{channel}/{syst}/pair_onZ/mass", pair.M(), totWeight, 60, 60., 120.)
            else:
                self.FillHist(f"{channel}/{syst}/pair_offZ/mass", pair.M(), totWeight, 200, 0., 200.)

        elif "2E1Mu" in channel:
            pair = electrons.at(0) + electrons.at(1)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), totWeight, 200, 0., 200.)
            # Delta R between leptons
            dR_ele1_mu = electrons.at(0).DeltaR(muons.at(0))
            dR_ele2_mu = electrons.at(1).DeltaR(muons.at(0))
            dR_ele1_ele2 = electrons.at(0).DeltaR(electrons.at(1))
            self.FillHist(f"{channel}/{syst}/dR_ele1_mu", dR_ele1_mu, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele2_mu", dR_ele2_mu, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele1_ele2", dR_ele1_ele2, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_min_ele_mu", min([dR_ele1_mu, dR_ele2_mu]), totWeight, 100, 0., 10.)

        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
            pair_lowM, pair_highM = (pair1, pair2) if pair1.M() < pair2.M() else (pair2, pair1)
            self.FillHist(f"{channel}/{syst}/pair_lowM/pt", pair_lowM.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair_lowM/eta", pair_lowM.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair_lowM/phi", pair_lowM.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair_lowM/mass", pair_lowM.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/pair_highM/pt", pair_highM.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair_highM/eta", pair_highM.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair_highM/phi", pair_highM.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair_highM/mass", pair_highM.M(), totWeight, 200, 0., 200.)

            if (60 < pair_lowM.M() and pair_lowM.M() < 120) or (60 < pair_highM.M() and pair_highM.M() < 120):
                self.FillHist(f"{channel}/{syst}/pair_lowM_onZ/mass", pair_lowM.M(), totWeight, 200, 0., 200.)
                self.FillHist(f"{channel}/{syst}/pair_highM_onZ/mass", pair_highM.M(), totWeight, 200, 0., 200.)
            else:
                self.FillHist(f"{channel}/{syst}/pair_lowM_offZ/mass", pair_lowM.M(), totWeight, 200, 0., 200.)
                self.FillHist(f"{channel}/{syst}/pair_highM_offZ/mass", pair_highM.M(), totWeight, 200, 0., 200.)

            # Delta R between leptons
            dR_pair_ss1_os = mu_ss1.DeltaR(mu_os)
            dR_pair_ss2_os = mu_ss2.DeltaR(mu_os)
            dR_pair_ss1_ss2 = mu_ss1.DeltaR(mu_ss2)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_os", dR_pair_ss1_os, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss2_os", dR_pair_ss2_os, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_ss2", dR_pair_ss1_ss2, totWeight, 100, 0., 10.)

        # Fill ZCands
        if "1E2Mu" in channel:
            ZCand = muons.at(0) + muons.at(1)
            nonprompt = electrons.at(0)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), totWeight, 64, -3.2, 3.2)

        if "2E1Mu" in channel:
            ZCand = electrons.at(0) + electrons.at(1)
            nonprompt = muons.at(0)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), totWeight, 64, -3.2, 3.2)

        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
            mZ = 91.2
            if abs(pair1.M() - mZ) < abs(pair2.M() - mZ):
                ZCand, nZCand = pair1, pair2
                nonprompt = mu_ss2
            else:
                ZCand, nZCand = pair2, pair1
                nonprompt = mu_ss1
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), totWeight, 64, -3.2, 3.2)

    def _fillCRObjects(self, channel: str,
                             recoObjects: dict,
                             weights: dict,
                             syst="Central"):
        """Fill CR-specific histograms (from CRPromptSelector)."""
        muons = recoObjects["tightMuons"]
        electrons = recoObjects["tightElectrons"]

        # Calculate total weight
        totWeight = weights["genWeight"] * weights["prefireWeight"] * weights["pileupWeight"]
        totWeight *= weights["muonRecoSF"] * weights["muonIDSF"] * weights["eleRecoSF"] * weights["eleIDSF"]
        totWeight *= weights["trigSF"] * weights["pileupIDSF"] * weights["btagSF"]
        totWeight *= weights["WZNjetsSF"]

        # Fill ZCand and conversion candidate for ZG channel
        if "ZG" in channel:
            ZCand = None
            if "1E2Mu" in channel:
                ZCand = electrons.at(0) + muons.at(0) + muons.at(1)
            if "3Mu" in channel:
                ZCand = muons.at(0) + muons.at(1) + muons.at(2)
            convLep = self.getConvLepton(muons, electrons)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/convLep/pt", convLep.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/convLep/eta", convLep.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/convLep/phi", convLep.Phi(), totWeight, 64, -3.2, 3.2)

        # Fill ZCand and nZCand for WZ channel
        if "WZ" in channel:
            ZCand = None
            if "1E2Mu" in channel:
                ZCand = muons.at(0) + muons.at(1)
                nZCand = muons.at(0) + muons.at(1)
                lep = electrons.at(0)
            if "3Mu" in channel:
                mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
                pair1, pair2 = (mu_ss1 + mu_os), (mu_ss2 + mu_os)
                if abs(pair1.M() - 91.2) < abs(pair2.M() - 91.2):
                    ZCand, nZCand = pair1, pair2
                    lep = mu_ss2
                else:
                    ZCand, nZCand = pair2, pair1
                    lep = mu_ss1
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), totWeight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/lep/pt", lep.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/lep/eta", lep.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/lep/phi", lep.Phi(), totWeight, 64, -3.2, 3.2)

    # ===================================================================
    # Tree Filling (for RunTreeMode)
    # ===================================================================

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
            # Create TTree using C++ NewTree (registers in treemap for WriteHist)
            thisTree = self.NewTree(f"Events_{syst}")

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

            self.trees[syst] = thisTree

    def fillTree(self, channel: str,
                       recoObjects: dict,
                       weights: dict,
                       syst="Central"):
        """Fill output tree with event information and ParticleNet scores."""
        if not self.RunTreeMode:
            return

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

        # ParticleNet scores (use scores from recoObjects)
        if syst == "Central":
            self.fold[syst][0] = recoObjects.get("fold", 0)
            scores = recoObjects.get("scores", {})
            for signal in self.signals:
                for cls in self.classNames:
                    score_key = f"{signal}_{cls}"
                    self.scores[syst][signal][cls][0] = scores.get(score_key, -999.)
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

    # Note: WriteHist is NOT overridden here.
    # Trees created via self.NewTree() are registered in C++ treemap
    # and written automatically by AnalyzerCore::WriteHist()


if __name__ == "__main__":
    module = PromptAnalyzer()
    module.SetTreeName("Events")
    module.LogEvery = 5000
    module.IsDATA = False
    module.MCSample = "TTLL_powheg"
    module.xsec = 98.03
    module.sumW = 3869737330.782318
    module.sumSign = 47713334.0
    module.SetEra("2023")
    module.SetPeriod("")
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst"])
    module.AddFile("tree_0.root")
    module.MaxEvent = max(1, int(module.fChain.GetEntries() / 1))
    module.SetOutfilePath("hists.root")
    module.Init()
    module.initializePyAnalyzer()
    module.Loop()
    module.WriteHist()
