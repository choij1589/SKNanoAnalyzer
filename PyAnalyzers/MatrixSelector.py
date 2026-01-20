import os
from ROOT import TString
from ROOT.VecOps import RVec
from ROOT import TriLeptonBase
from ROOT import JetTagging
from ROOT import Event, Lepton, Muon, Electron, Jet
from enum import IntEnum

from MLTools.helpers import loadMultiClassParticleNet, getGraphInput, getMultiClassScore

class CutStage(IntEnum):
    Initial = 0
    NoiseFilter = 1
    EventVetoMap = 2
    LeptonSelection = 3
    Trigger = 4
    KinematicCuts = 5
    JetRequirements = 6
    JetVetoMap = 7
    Final = 8

class MatrixSelector(TriLeptonBase):
    def __init__(self):
        super().__init__()

    def initializePyAnalyzer(self):
        self.initializeAnalyzer()

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
        self.signals = ["MHc160_MA85", "MHc130_MA90", "MHc130_MA100", "MHc100_MA95", "MHc115_MA87", "MHc145_MA92", "MHc160_MA98"]
        self.classNames = ["signal", "nonprompt", "diboson", "ttZ"]

        # Load ParticleNet models
        print(f"[MatrixSelector] Loading ParticleNet models for {self.channel}")
        self.models = loadMultiClassParticleNet(self.signals)
        print(f"[MatrixSelector] Loaded {len(self.models)} models")

    def fillCutflow(self, stage, channel, weight, syst):
        if not syst == "Central": return
        if weight is None: return
        cutIndex = int(stage)
        self.FillHist(f"{channel}/{syst}/cutflow", cutIndex, weight, 9, 0., 9.)

    def fillJetEtaPhi2D(self, jets, weight, stage):
        """Fill 2D jet eta-phi distribution for veto map validation."""
        for jet in jets:
            self.FillHist(f"JetEtaPhi/{stage}/eta_phi",
                          jet.Eta(), jet.Phi(), weight,
                          100, -5.0, 5.0,   # eta: 100 bins from -5 to 5
                          64, -3.2, 3.2)    # phi: 64 bins from -pi to pi

    def fillElectronScEtaPhi2D(self, electrons, weight, stage):
        """Fill 2D electron scEta-scPhi distribution for HEM veto validation."""
        for ele in electrons:
            self.FillHist(f"ElectronScEtaPhi/{stage}/scEta_scPhi",
                          ele.scEta(), ele.scPhi(), weight,
                          100, -2.5, 2.5,   # scEta: 100 bins from -2.5 to 2.5
                          64, -3.2, 3.2)    # scPhi: 64 bins from -pi to pi

    def executeEvent(self):
        ev = self.GetEvent()

        # Initial cutflow entry
        self.fillCutflow(CutStage.Initial, self.channel, 1.0, "Central")

        rawJets = self.GetAllJets()
        if not self.PassNoiseFilter(rawJets, ev): return
        self.fillCutflow(CutStage.NoiseFilter, self.channel, 1.0, "Central")

        rawMuons = self.GetAllMuons()
        if not (self.RunNoJetVeto or self.PassVetoMap(rawJets, rawMuons, "jetvetomap")): return
        self.fillCutflow(CutStage.EventVetoMap, self.channel, 1.0, "Central")

        # Fill jet eta-phi for events passing veto map (Run 3)
        if self.Run == 3:
            self.fillJetEtaPhi2D(rawJets, 1.0, "PassedEventVeto")

        rawElectrons = self.GetAllElectrons()

        recoObjects = self.defineObjects(ev, rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, recoObjects)
        if channel is None: return

        # Evaluate ParticleNet scores
        data, scores, fold = self.evalScore(
            recoObjects["looseMuons"],
            recoObjects["looseElectrons"],
            recoObjects["jets"],
            recoObjects["bjets"],
            recoObjects["METv"]
        )
        recoObjects["scores"] = scores

        weight = self.GetFakeWeight(recoObjects["looseMuons"], recoObjects["looseElectrons"], "Central")
        self.fillObjects(channel, recoObjects, weight, syst="Central")
    
    def defineObjects(self, ev, rawMuons, rawElectrons, rawJets):
        allMuons = RVec(Muon)(rawMuons)
        allElectrons = RVec(Electron)(rawElectrons)
        allJets = RVec(Jet)(rawJets)
        METv = ev.GetMETVector(Event.MET_Type.PUPPI)
        METv = self.ApplyTypeICorrection(METv, allJets, allElectrons, allMuons)

        # Sort objects in pt order
        allMuons = RVec(Muon)(sorted(allMuons, key=lambda x: x.Pt(), reverse=True))
        allElectrons = RVec(Electron)(sorted(allElectrons, key=lambda x: x.Pt(), reverse=True))
        allJets = RVec(Jet)(sorted(allJets, key=lambda x: x.Pt(), reverse=True))

        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        looseMuons = self.SelectMuons(vetoMuons, self.MuonIDs.GetID("loose"), 10., 2.4)
        tightMuons = self.SelectMuons(looseMuons, self.MuonIDs.GetID("tight"), 10., 2.4)

        # HEM veto: only apply to loose/tight electrons, not veto (to properly reject events with extra leptons)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs.GetID("loose"), 10., 2.5)
        applyHEMVeto = (self.DataEra == "2018") and not self.RunNoHEMVeto

        # Fill electron scEta-scPhi BEFORE HEM veto (2018 only)
        if applyHEMVeto:
            self.fillElectronScEtaPhi2D(vetoElectrons, 1.0, "BeforeHEMVeto")

        looseElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs.GetID("loose"), 15., 2.5, applyHEMVeto)
        tightElectrons = self.SelectElectrons(looseElectrons, self.ElectronIDs.GetID("tight"), 15., 2.5, applyHEMVeto)

        # Fill electron scEta-scPhi AFTER HEM veto (2018 only)
        if applyHEMVeto:
            self.fillElectronScEtaPhi2D(looseElectrons, 1.0, "AfterHEMVeto")

        max_jeteta = 2.4 if self.DataEra.Contains("2016") else 2.5
        jets_passtight = self.SelectJets(allJets, "tight", 20., max_jeteta)
        jets_vetoLep = self.JetsVetoLeptonInside(jets_passtight, vetoElectrons, vetoMuons, 0.4)

        jets = RVec(Jet)()
        bjets = RVec(Jet)()
        jets_beforeVeto = RVec(Jet)()  # Track jets before veto map (Run 2 only)
        tagger = JetTagging.JetFlavTagger.DeepJet
        wp = self.myCorr.GetBTaggingWP(tagger, JetTagging.JetFlavTaggerWP.Medium)

        for j in jets_vetoLep:
            if self.Run == 2:
                if not j.PassID("loosePuId"): continue
                jets_beforeVeto.emplace_back(j)  # Track before veto map
                if not (self.RunNoJetVeto or self.PassVetoMap(j, allMuons, "jetvetomap")): continue
            jets.emplace_back(j)

            if j.GetBTaggerResult(tagger) > wp:
                bjets.emplace_back(j)

        # Fill jet eta-phi for Run 2
        if self.Run == 2:
            self.fillJetEtaPhi2D(jets_beforeVeto, 1.0, "BeforeJetVeto")
            self.fillJetEtaPhi2D(jets, 1.0, "AfterJetVeto")


        return {"vetoMuons": vetoMuons,
                "looseMuons": looseMuons,
                "tightMuons": tightMuons,
                "vetoElectrons": vetoElectrons,
                "looseElectrons": looseElectrons,
                "tightElectrons": tightElectrons,
                "jets": jets,
                "bjets": bjets,
                "METv": METv}
    
    def selectEvent(self, ev: Event,
                          recoObjects: dict,
                          weight: float = 1.0) -> str:
        vetoMuons = recoObjects["vetoMuons"]
        looseMuons = recoObjects["looseMuons"]
        tightMuons = recoObjects["tightMuons"]
        vetoElectrons = recoObjects["vetoElectrons"]
        looseElectrons = recoObjects["looseElectrons"]
        tightElectrons = recoObjects["tightElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]

        is1E2Mu = (looseElectrons.size() == 1 and vetoElectrons.size() == 1 and \
                   looseMuons.size() == 2 and vetoMuons.size() == 2)
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is2E1Mu = (looseElectrons.size() == 2 and vetoElectrons.size() == 2 and \
                   looseMuons.size() == 1 and vetoMuons.size() == 1)

        #### Not all leptons tight
        if self.Run1E2Mu:
            if not is1E2Mu: return
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return
            self.fillCutflow(CutStage.LeptonSelection, self.channel, weight, "Central")

        if self.Run3Mu:
            if not is3Mu: return
            if tightMuons.size() == looseMuons.size(): return
            self.fillCutflow(CutStage.LeptonSelection, self.channel, weight, "Central")

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers): return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2 = looseMuons.at(0), looseMuons.at(1)
            ele = looseElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return

            # Scale loose leptons after trigger cuts and update recoObjects
            looseMuons = self.GetPTCorrScaledMuons(looseMuons)
            looseElectrons = self.GetPTCorrScaledElectrons(looseElectrons)
            recoObjects["looseMuons"] = looseMuons
            recoObjects["looseElectrons"] = looseElectrons

            if not mu1.Charge()+mu2.Charge() == 0: return

            pair = mu1 + mu2
            if not pair.M() > 12.: return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if not jets.size() >= 2: return
            self.fillCutflow(CutStage.JetRequirements, self.channel, weight, "Central")

            if bjets.size() == 0:
                isOnZ = abs(pair.M() - 91.2) < 10.
                if isOnZ:
                    self.fillCutflow(CutStage.Final, "ZFake1E2Mu", weight, "Central")
                    return "ZFake1E2Mu"
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
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            mu1, mu2, mu3 = tuple(looseMuons)
            if not mu1.Pt() > 20.: return
            if not mu2.Pt() > 10.: return
            if not mu3.Pt() > 10.: return

            # Scale loose muons after trigger cuts and update recoObjects
            looseMuons = self.GetPTCorrScaledMuons(looseMuons)
            recoObjects["looseMuons"] = looseMuons

            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(looseMuons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            if not pair1.M() > 12.: return
            if not pair2.M() > 12.: return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if not jets.size() >= 2: return
            self.fillCutflow(CutStage.JetRequirements, self.channel, weight, "Central")

            if bjets.size() == 0:
                isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
                if isOnZ:
                    self.fillCutflow(CutStage.Final, "ZFake3Mu", weight, "Central")
                    return "ZFake3Mu"
                return
            else:
                self.fillCutflow(CutStage.Final, "SR3Mu", weight, "Central")
                return "SR3Mu"

        ## 2E1Mu sideband (TTZ control region)
        ## 1. pass EMuTriggers
        ## 2. Exact 2 loose electrons and 1 loose muon, no additional leptons
        ## 3. NOT all leptons tight (matrix method requirement)
        ## 4. Exist OS electron pair, 60 < M(ee) < 120 GeV
        ## 5. At least two jets, at least one bjet
        if self.Run2E1Mu:
            if not is2E1Mu: return
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return
            self.fillCutflow(CutStage.LeptonSelection, self.channel, weight, "Central")

            if not ev.PassTrigger(self.EMuTriggers): return
            self.fillCutflow(CutStage.Trigger, self.channel, weight, "Central")

            el1, el2 = looseElectrons.at(0), looseElectrons.at(1)
            mu = looseMuons.at(0)
            passLeadMu = mu.Pt() > 25. and el1.Pt() > 15. and el2.Pt() > 15.
            passLeadEl = (el1.Pt() > 25. or el2.Pt() > 25.) and mu.Pt() > 10.
            passSafeCut = passLeadMu or passLeadEl
            if not passSafeCut: return

            # Scale loose leptons after trigger cuts and update recoObjects
            looseMuons = self.GetPTCorrScaledMuons(looseMuons)
            looseElectrons = self.GetPTCorrScaledElectrons(looseElectrons)
            recoObjects["looseMuons"] = looseMuons
            recoObjects["looseElectrons"] = looseElectrons

            if not (el1.Charge() + el2.Charge() == 0): return

            el1, el2 = looseElectrons.at(0), looseElectrons.at(1)
            pair = el1 + el2
            if not 60. < pair.M() < 120.: return
            self.fillCutflow(CutStage.KinematicCuts, self.channel, weight, "Central")

            if not jets.size() >= 2: return
            self.fillCutflow(CutStage.JetRequirements, self.channel, weight, "Central")

            if not bjets.size() >= 1: return
            self.fillCutflow(CutStage.Final, "TTZ2E1Mu", weight, "Central")
            return "TTZ2E1Mu"
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

    def fillObjects(self, channel: str,
                          recoObjects: dict,
                          weight: float,
                          syst: str = "Central"):
        electrons = recoObjects["looseElectrons"]
        muons = recoObjects["looseMuons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]
        
        ## fill base observables
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/energy", mu.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/px", mu.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/py", mu.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pz", mu.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/charge", mu.Charge(), weight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/muons/{idx}/min_dR_bjets", min([mu.DeltaR(bjet) for bjet in bjets]), weight, 100, 0., 10.)
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/energy", ele.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/px", ele.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/py", ele.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pz", ele.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/charge", ele.Charge(), weight, 3, -1, 2)
            if bjets.size() > 0:
                self.FillHist(f"{channel}/{syst}/electrons/{idx}/min_dR_bjets", min([ele.DeltaR(bjet) for bjet in bjets]), weight, 100, 0., 10.)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/energy", jet.E(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/px", jet.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/py", jet.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pz", jet.Pz(), weight, 500, -250., 250.)
        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/energy", bjet.E(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/px", bjet.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/py", bjet.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pz", bjet.Pz(), weight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/{syst}/METv/px", METv.Px(), weight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/py", METv.Py(), weight, 500, -250., 250.)

        # Fill discrimination variable
        if "1E2Mu" in channel:
            pair = muons.at(0) + muons.at(1)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 200, 0., 200.)
            ## Delta R between leptons
            dR_ele_mu1 = electrons.at(0).DeltaR(muons.at(0))
            dR_ele_mu2 = electrons.at(0).DeltaR(muons.at(1))
            dR_mu1_mu2 = muons.at(0).DeltaR(muons.at(1))
            self.FillHist(f"{channel}/{syst}/dR_ele_mu1", dR_ele_mu1, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele_mu2", dR_ele_mu2, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_mu1_mu2", dR_mu1_mu2, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_min_ele_mu", min([dR_ele_mu1, dR_ele_mu2]), weight, 100, 0., 10.)

            if 60 < pair.M() and pair.M() < 120:
                self.FillHist(f"{channel}/{syst}/pair_onZ/mass", pair.M(), weight, 60, 60., 120.)
            else:
                self.FillHist(f"{channel}/{syst}/pair_offZ/mass", pair.M(), weight, 200, 0., 200.)
        elif "2E1Mu" in channel:
            pair = electrons.at(0) + electrons.at(1)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 200, 0., 200.)
            ## Delta R between leptons
            dR_ele1_mu = electrons.at(0).DeltaR(muons.at(0))
            dR_ele2_mu = electrons.at(1).DeltaR(muons.at(0))
            dR_ele1_ele2 = electrons.at(0).DeltaR(electrons.at(1))
            self.FillHist(f"{channel}/{syst}/dR_ele1_mu", dR_ele1_mu, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele2_mu", dR_ele2_mu, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_ele1_ele2", dR_ele1_ele2, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_min_ele_mu", min([dR_ele1_mu, dR_ele2_mu]), weight, 100, 0., 10.)
        elif "3Mu" in channel:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            pair_lowM, pair_highM = (pair1, pair2) if pair1.M() < pair2.M() else (pair2, pair1)
            self.FillHist(f"{channel}/{syst}/pair_lowM/pt", pair_lowM.Pt(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair_lowM/eta", pair_lowM.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair_lowM/phi", pair_lowM.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair_lowM/mass", pair_lowM.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/pair_highM/pt", pair_highM.Pt(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/pair_highM/eta", pair_highM.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair_highM/phi", pair_highM.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair_highM/mass", pair_highM.M(), weight, 200, 0., 200.)

            dR_pair_ss1_os = mu_ss1.DeltaR(mu_os)
            dR_pair_ss2_os = mu_ss2.DeltaR(mu_os)
            dR_pair_ss1_ss2 = mu_ss1.DeltaR(mu_ss2)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_os", dR_pair_ss1_os, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss2_os", dR_pair_ss2_os, weight, 100, 0., 10.)
            self.FillHist(f"{channel}/{syst}/dR_pair_ss1_ss2",  dR_pair_ss1_ss2, weight, 100, 0., 10.)

            if (60 < pair_lowM.M() and pair_lowM.M() < 120) or (60 < pair_highM.M() and pair_highM.M() < 120):
                self.FillHist(f"{channel}/{syst}/pair_lowM_onZ/mass", pair_lowM.M(), weight, 200, 0., 200.)
                self.FillHist(f"{channel}/{syst}/pair_highM_onZ/mass", pair_highM.M(), weight, 200, 0., 200.)
            else:
                self.FillHist(f"{channel}/{syst}/pair_lowM_offZ/mass", pair_lowM.M(), weight, 200, 0., 200.)
                self.FillHist(f"{channel}/{syst}/pair_highM_offZ/mass", pair_highM.M(), weight, 200, 0., 200.)

        # Fill ZCands
        if "1E2Mu" in channel:
            ZCand = muons.at(0) + muons.at(1)
            nonprompt = electrons.at(0)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
        elif "2E1Mu" in channel:
            ZCand = electrons.at(0) + electrons.at(1)
            nonprompt = muons.at(0)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
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
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 500, 0., 500.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)

        # Fill ParticleNet scores
        if "scores" in recoObjects:
            scores = recoObjects["scores"]
            for signal in self.signals:
                # Fill multi-class ParticleNet scores
                score_signal    = scores[f"{signal}_signal"]
                score_nonprompt = scores[f"{signal}_nonprompt"]
                score_diboson   = scores[f"{signal}_diboson"]
                score_ttZ       = scores[f"{signal}_ttZ"]
                self.FillHist(f"{channel}/{syst}/{signal}/score_signal", score_signal, weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_nonprompt", score_nonprompt, weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_diboson", score_diboson, weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/score_ttZ", score_ttZ, weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/LR_nonprompt", score_signal/(score_signal+score_nonprompt), weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/LR_diboson", score_signal/(score_signal+score_diboson), weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/LR_ttZ", score_signal/(score_signal+score_ttZ), weight, 100, 0., 1.)
                self.FillHist(f"{channel}/{syst}/{signal}/LR_totalBkg", score_signal/(score_signal+score_nonprompt+score_diboson+score_ttZ), weight, 100, 0., 1.)

if __name__ == "__main__":
    module = MatrixSelector()
    module.SetTreeName("Events")
    module.LogEvery = 5000
    module.IsDATA = True
    module.DataStream = "MuonEG"
    module.SetEra("2017")
    module.SetPeriod("C")
    module.Userflags = RVec(TString)(["Run1E2Mu", "RunSyst"])
    module.AddFile("tree_0.root")
    module.MaxEvent = max(1, int(module.fChain.GetEntries()/1))
    module.SetOutfilePath("hists_0.root")
    module.Init()
    module.initializePyAnalyzer()
    module.Loop()
    module.WriteHist()
