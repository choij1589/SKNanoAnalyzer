import os
from ROOT import TString
from ROOT.VecOps import RVec
from ROOT import TriLeptonBase
from ROOT import MyCorrection; myVar = MyCorrection.variation
from ROOT import SystematicHelper
from ROOT import JetTagging
from ROOT import Event, Lepton, Muon, Electron, Jet
from ROOT import Gen, GenJet
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

class PromptSelector(TriLeptonBase):
    def __init__(self):
        super().__init__()
        self.systHelper = None
    
    def fillCutflow(self, stage, channel, weight, syst):
        if not syst == "Central": return
        if weight is None: return
        channelStr = channel if channel else "ALL"
        cutIndex = int(stage)
        self.FillHist(f"{channelStr}/{syst}/cutflow", cutIndex, weight, 10, 0., 10.)
        
    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
        
        # Link flags
        if not (self.Run1E2Mu or self.Run3Mu):
            raise ValueError("Run1E2Mu or Run3Mu must be set")
        if self.Run1E2Mu and self.Run3Mu:
            raise ValueError("Run1E2Mu and Run3Mu cannot be set at the same time")
        
        # Determine channel
        self.channel = "Run1E2Mu" if self.Run1E2Mu else "Run3Mu"

        # ParticleNet configuration
        self.signals = ["MHc160_MA85", "MHc130_MA90", "MHc130_MA100", "MHc100_MA95", "MHc115_MA87", "MHc145_MA92", "MHc160_MA98"]
        self.classNames = ["signal", "nonprompt", "diboson", "ttZ"]

        # Load ParticleNet models
        print(f"[PromptSelector] Loading ParticleNet models for {self.channel}")
        self.models = loadMultiClassParticleNet(self.signals)
        print(f"[PromptSelector] Loaded {len(self.models)} models")

        # Systematics
        if self.IsDATA:
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml", self.DataStream, self.DataEra)
        elif not self.RunSyst: # MC without RunSyst
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/noSyst.yaml", self.MCSample, self.DataEra)
        else: # MC and RunSyst
            self.systHelper = SystematicHelper(f"{os.getenv('SKNANO_HOME')}/AnalyzerTools/TriLeptonSystematics.yaml", self.MCSample, self.DataEra)
        
    def executeEvent(self):
        ev = self.GetEvent()
        
        # Initial cutflow entry
        initialWeight = 1.0 if self.IsDATA else self.MCweight() * ev.GetTriggerLumi("Full")
        self.fillCutflow(CutStage.Initial, None, initialWeight, "Central")
        
        rawJets = self.GetAllJets()
        if not self.PassNoiseFilter(rawJets, ev): return
        self.fillCutflow(CutStage.NoiseFilter, None, initialWeight, "Central")
        
        rawMuons = self.GetAllMuons()
        if not (self.RunNoVetoMap or self.PassVetoMap(rawJets, rawMuons, "jetvetomap")): return
        self.fillCutflow(CutStage.VetoMap, None, initialWeight, "Central")
        
        rawElectrons = self.GetAllElectrons()
        truth = self.GetAllGens()
        genJets = self.GetAllGenJets()

        # Cache base object copies once per event
        self._cached_objects = {}
        
        def processEvent(syst, apply_weight_variation=False):
            recoObjects = self.defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, syst)
            channel = self.selectEvent(ev, recoObjects, truth, syst, initialWeight if syst == "Central" else None)
            if channel is None: return

            # Evaluate ParticleNet scores
            data, scores, fold = self.evalScore(
                recoObjects["tightMuons"],
                recoObjects["tightElectrons"],
                recoObjects["jets"],
                recoObjects["bjets"],
                recoObjects["METv"]
            )
            recoObjects["scores"] = scores

            if apply_weight_variation:
                assert syst == "Central", "Only Central weight variation is allowed"
                weights = self.getWeights(ev, recoObjects, genJets, "Central")
                self.fillObjects(channel, recoObjects, weights, "Central")

                ## Get weight-only systematics
                weightOnlySysts = list(self.systHelper.getWeightOnlySystematics())
                if self.RunSyst and "1E2Mu" in channel:
                    weightOnlySysts.remove("DblMuTrigSF")
                if self.RunSyst and "3Mu" in channel:
                    weightOnlySysts.remove("ElectronIDSF")
                    weightOnlySysts.remove("EMuTrigSF")

                for syst in weightOnlySysts:
                    weights_up = self.getWeights(ev, recoObjects, genJets, f"{syst}_Up")
                    weights_down = self.getWeights(ev, recoObjects, genJets, f"{syst}_Down")
                    self.fillObjects(channel, recoObjects, weights_up, f"{syst}_Up")
                    self.fillObjects(channel, recoObjects, weights_down, f"{syst}_Down")
            else:
                weights = self.getWeights(ev, recoObjects, genJets, syst)
                self.fillObjects(channel, recoObjects, weights, syst)
        
        processEvent("Central", apply_weight_variation=True)
        
        # Process systematics requiring evtLoopAgain
        for syst in self.systHelper:
            systName = syst.iter_name
            
            # Skip Central (already processed) and weight-only systematics
            if systName == "Central": continue
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
            variation = syst.split("_")[1].lower() # up or down
            allElectrons = self.ScaleElectrons(ev, allElectrons, variation)
        if "ElectronRes" in syst:
            variation = syst.split("_")[1].lower() # up or down
            allElectrons = self.SmearElectrons(allElectrons, variation)
        if "MuonEn" in syst:
            variation = syst.split("_")[1].lower() # up or down
            allMuons = self.ScaleMuons(allMuons, variation)
        if "JetEn" in syst:
            variation = syst.split("_")[1].lower() # up or down
            allJets = self.ScaleJets(allJets, variation, "total")
        if "JetRes" in syst:
            variation = syst.split("_")[1].lower() # up or down
            allJets = self.SmearJets(allJets, genJets, variation)
        MET_var = Event.MET_Syst.CENTRAL
        if "UnclusteredEn" in syst:
            if syst.split("_")[1].lower() == "up":
                MET_var = Event.MET_Syst.UE_UP
            else:
                MET_var = Event.MET_Syst.UE_DOWN
        METv = ev.GetMETVector(Event.MET_Type.PUPPI, MET_var)
        METv = self.ApplyTypeICorrection(METv, allJets, allElectrons, allMuons)
        
        # sort objects in pt order
        allMuons = RVec(Muon)(sorted(allMuons, key=lambda x: x.Pt(), reverse=True))
        allElectrons = RVec(Electron)(sorted(allElectrons, key=lambda x: x.Pt(), reverse=True))
        allJets = RVec(Jet)(sorted(allJets, key=lambda x: x.Pt(), reverse=True))
        
        # select objects
        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        tightMuons = self.SelectMuons(vetoMuons, self.MuonIDs.GetID("tight"), 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs.GetID("loose"), 10., 2.5)
        tightElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs.GetID("tight"), 15., 2.5)
        
        max_jeteta = 2.4 if self.DataEra.Contains("2016") else 2.5
        jets_selected = self.SelectJets(allJets, "tight", 20., max_jeteta)
        jets_vetoLep = self.JetsVetoLeptonInside(jets_selected, vetoElectrons, vetoMuons, 0.4)
        
        # Single-pass filtering for better performance
        jets = RVec(Jet)()
        bjets = RVec(Jet)()
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = self.myCorr.GetBTaggingWP(tagger, JetTagging.JetFlavTaggerWP.Medium)
        
        for j in jets_vetoLep:
            # Apply Run 2 specific filters
            if self.Run == 2:
                if not j.PassID("loosePuId"): continue
                if not (self.RunNoVetoMap or self.PassVetoMap(j, allMuons, "jetvetomap")): continue
            jets.emplace_back(j)
            
            # b-tagging
            if j.GetBTaggerResult(tagger) > wp:
                bjets.emplace_back(j)

        return {"vetoMuons": vetoMuons,
                "tightMuons": tightMuons,
                "vetoElectrons": vetoElectrons,
                "tightElectrons": tightElectrons,
                "jets_vetoLep": jets_vetoLep,
                "jets": jets,
                "bjets": bjets,
                "METv": METv}

    def selectEvent(self, ev: Event, 
                          recoObjects: dict,
                          truth: RVec[Gen], 
                          syst: str = "Central",
                          weight: float = None) -> str:
        vetoMuons = recoObjects["vetoMuons"]
        tightMuons = recoObjects["tightMuons"]
        vetoElectrons = recoObjects["vetoElectrons"]
        tightElectrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        is1E2Mu = (tightElectrons.size() == 1 and vetoElectrons.size() == 1 and \
                   tightMuons.size() == 2 and vetoMuons.size() == 2)
        is3Mu = (tightMuons.size() == 3 and vetoMuons.size() == 3 and \
                 tightElectrons.size() == 0 and vetoElectrons.size() == 0)
        
        if self.Run1E2Mu and not is1E2Mu: return
        if self.Run3Mu and not is3Mu: return

        # Record lepton selection cutflow 
        if is1E2Mu:
            self.fillCutflow(CutStage.LeptonSelection, "1E2Mu", weight, "Central")
        elif is3Mu:
            self.fillCutflow(CutStage.LeptonSelection, "3Mu", weight, "Central")

        # for conversion samples
        if self.MCSample.Contains("DYJets") or self.MCSample.Contains("ZGToLLG") or self.MCSample.Contains("DYGTo2LG"):
            # at least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = RVec(Muon)()
            fakeMuons = RVec(Muon)()
            convElectrons = RVec(Electron)()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
                if self.GetLeptonType(mu, truth) in [-1, -2, -3, -4]: fakeMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
            # remove hadronic contribution
            if self.Run1E2Mu:
                if not fakeMuons.size() == 0: return
                if not (convElectrons.size()+convMuons.size()) > 0: return
                self.fillCutflow(CutStage.ConversionFilter, "1E2Mu", weight, "Central")
            if self.Run3Mu:
                if not fakeMuons.size() == 0: return
                if not convMuons.size() > 0: return
                self.fillCutflow(CutStage.ConversionFilter, "3Mu", weight, "Central")
        # Patching sample
        leptons = RVec(Lepton)()
        for mu in tightMuons: leptons.emplace_back(mu)
        for ele in tightElectrons: leptons.emplace_back(ele)
        region = "LowPT" if any(l.Pt() < 15. for l in leptons) else "HighPT"
        if self.MCSample.Contains("DYJets") and not region == "LowPT": return
        if (self.MCSample.Contains("ZGToLLG") or self.MCSample.Contains("DYGTo2LG")) and not region == "HighPT": return
        
        # Sample patching cutflow
        if syst == "Central" and weight is not None:
            if self.Run1E2Mu:
                self.fillCutflow(CutStage.SamplePatching, "1E2Mu", weight, "Central")
            elif self.Run3Mu:
                self.fillCutflow(CutStage.SamplePatching, "3Mu", weight, "Central")

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers): return
            if syst == "Central" and weight is not None:
                self.fillCutflow(CutStage.Trigger, "1E2Mu", weight, "Central")
                
            mu1, mu2 = tightMuons.at(0), tightMuons.at(1)
            ele = tightElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return
            if not mu1.Charge()+mu2.Charge() == 0: return
            pair = mu1 + mu2
            if not pair.M() > 12.: return
            self.fillCutflow(CutStage.KinematicCuts, "1E2Mu", weight, "Central")
            
            if jets.size() >= 2:
                self.fillCutflow(CutStage.JetRequirements, "1E2Mu", weight, "Central")
                if bjets.size() == 0:
                    isOnZ = abs(pair.M() - 91.2) < 10.
                    if isOnZ:
                        self.fillCutflow(CutStage.Final, "ZFake1E2Mu", weight, "Central")
                        return "ZFake1E2Mu"
                    else:
                        return
                else:
                    self.fillCutflow(CutStage.Final, "SR1E2Mu", weight, "Central")
                    return "SR1E2Mu"
        
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets
        if self.Run3Mu:
            if not ev.PassTrigger(self.DblMuTriggers): return
            self.fillCutflow(CutStage.Trigger, "3Mu", weight, "Central")
                
            mu1, mu2, mu3 = tuple(tightMuons)
            if not mu1.Pt() > 20.: return
            if not mu2.Pt() > 10.: return
            if not mu3.Pt() > 10.: return
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            if not pair1.M() > 12.: return
            if not pair2.M() > 12.: return
            self.fillCutflow(CutStage.KinematicCuts, "3Mu", weight, "Central")
                
            if not jets.size() >= 2: return
            if syst == "Central" and weight is not None:
                self.fillCutflow(CutStage.JetRequirements, "3Mu", weight, "Central")
            if bjets.size() == 0:
                isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
                if isOnZ: 
                    self.fillCutflow(CutStage.Final, "ZFake3Mu", weight, "Central")
                    return "ZFake3Mu"
                return
            else:
                self.fillCutflow(CutStage.Final, "SR3Mu", weight, "Central")
                return "SR3Mu"
        return 

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

    def evalScore(self, muons, electrons, jets, bjets, METv):
        """Evaluate ParticleNet scores for all signal mass points."""
        scores = {}
        data, fold = getGraphInput(muons, electrons, jets, bjets, METv, str(self.DataEra))

        for sig in self.signals:
            if sig not in self.models.keys():
                print(f"[WARNING] Model {model_key} not found!")
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
        
        genWeight = self.MCweight()*ev.GetTriggerLumi("Full")

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
        else:
            pass
        
        btagSF = 1.
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = JetTagging.JetFlavTaggerWP.Medium
        method = JetTagging.JetTaggingSFMethod.mujets
        source = "central"
        if "HFcorr" in syst: source = "hf_corr"
        if "HFuncorr" in syst: source = "hf_uncorr"
        if "LFcorr" in syst: source = "lf_corr"
        if "LFuncorr" in syst: source = "lf_uncorr"
        btagSF = self.myCorr.GetBTaggingReweightMethod1a(recoObjects["jets"], tagger, wp, method, var, source)
        
        WZNjetsSF = 1.
        if (self.Run == 3 and (self.MCSample.Contains("WZTo3LNu") or self.MCSample.Contains("ZZTo4L"))) and (not self.RunNoWZSF):
            njets = float(recoObjects["jets"].size())
            WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "Central")
            if "WZNjetsSF_prompt" in syst:
                if "Up" in syst:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "prompt_up")
                elif "Down" in syst:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "prompt_down")
            if "WZNjetsSF_nonprompt" in syst:
                if "Up" in syst:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "nonprompt_up")
                else:
                    WZNjetsSF = self.myCorr.GetWZNjetsSF(njets, "nonprompt_down")
        
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
    
    def fillObjects(self, channel: str,
                          recoObjects: dict,
                          weights: dict,
                          syst="Central"):
        muons = recoObjects["tightMuons"]
        electrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        ## fill weights
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
        totWeight = genWeight*prefireWeight*pileupWeight
        totWeight *= muonRecoSF*muonIDSF*eleRecoSF*eleIDSF
        totWeight *= trigSF*pileupIDSF*btagSF
        totWeight *= WZNjetsSF

        ## fill weights
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

        ## fill base observables
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), totWeight, 10, 0., 1.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/energy", mu.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/px", mu.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/py", mu.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pz", mu.Pz(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/charge", mu.Charge(), totWeight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/muons/{idx}/min_dR_bjets", min([mu.DeltaR(bjet) for bjet in bjets]), totWeight, 100, 0., 10.)
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), totWeight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/energy", ele.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/px", ele.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/py", ele.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pz", ele.Pz(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/charge", ele.Charge(), totWeight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/electrons/{idx}/min_dR_bjets", min([ele.DeltaR(bjet) for bjet in bjets]), totWeight, 100, 0., 10.)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), totWeight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/energy", jet.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/px", jet.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/py", jet.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pz", jet.Pz(), totWeight, 500, -250., 250.)

        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), totWeight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/energy", bjet.E(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/px", bjet.Px(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/py", bjet.Py(), totWeight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pz", bjet.Pz(), totWeight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), totWeight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), totWeight, 15, 0., 15.)
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), totWeight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), totWeight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/{syst}/METv/energy", METv.E(), totWeight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/px", METv.Px(), totWeight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/py", METv.Py(), totWeight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/pz", METv.Pz(), totWeight, 500, -250., 250.)

        if "1E2Mu" in channel:
            pair = muons.at(0) + muons.at(1)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), totWeight, 200, 0., 200.)
            ## Delta R between leptons
            dR_ele_mu1 = electrons.at(0).DeltaR(muons.at(0))
            dR_ele_mu2 = electrons.at(0).DeltaR(muons.at(1))
            dR_mu1_mu2 = muons.at(0).DeltaR(muons.at(1))
            self.FillHist(f"{channel}/{syst}/dR_ele_mu1", dR_ele_mu1, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele_mu2", dR_ele_mu2, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_mu1_mu2", dR_mu1_mu2, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_min_ele_mu", min([dR_ele_mu1, dR_ele_mu2]), totWeight, 100, 0., 10.)

            # Within Z mass window?
            if 60 < pair.M() and pair.M() < 120:
                self.FillHist(f"{channel}/{syst}/pair_onZ/mass", pair.M(), totWeight, 60, 60., 120.);
            else:
                self.FillHist(f"{channel}/{syst}/pair_offZ/mass", pair.M(), totWeight, 200, 0., 200.);
        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            pair_lowM, pair_highM = (pair1, pair2) if pair1.M() < pair2.M() else (pair2, pair1)
            self.FillHist(f"{channel}/{syst}/pair_lowM/pt", pair_lowM.Pt(), totWeight, 300, 0., 300.);
            self.FillHist(f"{channel}/{syst}/pair_lowM/eta", pair_lowM.Eta(), totWeight, 100, -5., 5.);
            self.FillHist(f"{channel}/{syst}/pair_lowM/phi", pair_lowM.Phi(), totWeight, 64, -3.2, 3.2);
            self.FillHist(f"{channel}/{syst}/pair_lowM/mass", pair_lowM.M(), totWeight, 200, 0., 200.);
            self.FillHist(f"{channel}/{syst}/pair_highM/pt", pair_highM.Pt(), totWeight, 300, 0., 300.);
            self.FillHist(f"{channel}/{syst}/pair_highM/eta", pair_highM.Eta(), totWeight, 100, -5., 5.);
            self.FillHist(f"{channel}/{syst}/pair_highM/phi", pair_highM.Phi(), totWeight, 64, -3.2, 3.2);
            self.FillHist(f"{channel}/{syst}/pair_highM/mass", pair_highM.M(), totWeight, 200, 0., 200.);

            if (60 < pair_lowM.M() and pair_lowM.M() < 120) or (60 < pair_highM.M() and pair_highM.M() < 120):
                self.FillHist(f"{channel}/{syst}/pair_lowM_onZ/mass", pair_lowM.M(), totWeight, 200, 0., 200.);
                self.FillHist(f"{channel}/{syst}/pair_highM_onZ/mass", pair_highM.M(), totWeight, 200, 0., 200.);
            else:
                self.FillHist(f"{channel}/{syst}/pair_lowM_offZ/mass", pair_lowM.M(), totWeight, 200, 0., 200.);
                self.FillHist(f"{channel}/{syst}/pair_highM_offZ/mass", pair_highM.M(), totWeight, 200, 0., 200.);

            ## Delta R between leptons
            dR_pair_ss1_os = mu_ss1.DeltaR(mu_os)
            dR_pair_ss2_os = mu_ss2.DeltaR(mu_os)
            dR_pair_ss1_ss2 = mu_ss1.DeltaR(mu_ss2)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_os", dR_pair_ss1_os, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss2_os", dR_pair_ss2_os, totWeight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_ss2",  dR_pair_ss1_ss2, totWeight, 100, 0., 10.)
        
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
        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            mZ = 91.2
            if abs(pair1.M() - mZ) < abs(pair2.M() - mZ): 
                ZCand, nZCand = pair1, pair2
                nonprompt = mu_ss2
            else:                                         
                ZCand, nZCand = pair2, pair1
                nonprompt = mu_ss1
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), totWeight, 64, -3.2, 3.2)

        # Fill ParticleNet scores
        if "scores" in recoObjects:
            scores = recoObjects["scores"]
            for signal in self.signals:
                # Fill multi-class ParticleNet scores
                score_signal    = scores[f"{signal}_signal"]
                score_nonprompt = scores[f"{signal}_nonprompt"]
                score_diboson   = scores[f"{signal}_diboson"]
                score_ttZ       = scores[f"{signal}_ttZ"]
                self.FillHist(f"{channel}/{syst}/{signal}/score_signal", score_signal, totWeight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_nonprompt", score_nonprompt, totWeight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_diboson", score_diboson, totWeight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_ttZ", score_ttZ, totWeight, 100, 0., 1.)

if __name__ == "__main__":
    module = PromptSelector()
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
    module.MaxEvent = max(1, int(module.fChain.GetEntries()/1))
    module.SetOutfilePath("hists.root")
    module.Init()
    module.initializePyAnalyzer()
    module.Loop()
    module.WriteHist()
