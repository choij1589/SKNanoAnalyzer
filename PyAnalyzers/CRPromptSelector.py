import os
from ROOT import TString
from ROOT.VecOps import RVec
from PromptSelector import PromptSelector, CutStage
from ROOT import MyCorrection; myVar = MyCorrection.variation
from ROOT import Event, Lepton, Muon, Electron
from ROOT import Gen

class CRPromptSelector(PromptSelector):
    def __init__(self):
        super().__init__()
    
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
                if not convElectrons.size() == 1: return
                self.fillCutflow(CutStage.ConversionFilter, "1E2Mu", weight, "Central")
            if self.Run3Mu:
                if not fakeMuons.size() == 0: return
                if not convMuons.size() == 1: return
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

        # 1E2Mu ZGamma
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. |M(mumue) - 91.2| < 10
        ## 5. No b-jet
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
            if not (pair.M() > 12.): return
            self.fillCutflow(CutStage.KinematicCuts, "1E2Mu", weight, "Central")
            
            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair.M() - 91.2) < 10.): return
                if not bjets.size() == 0: return
                self.fillCutflow(CutStage.Final, "WZ1E2Mu", weight, "Central")
                return "WZ1E2Mu"
            else:
                if not (abs(pair.M() - 91.2) > 10.): return
                if not abs((mu1+mu2+ele).M() - 91.2) < 10.: return
                if not bjets.size() == 0: return
                self.fillCutflow(CutStage.Final, "ZG1E2Mu", weight, "Central")
                return "ZG1E2Mu"
        
        ## 3Mu ZGamma
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. |M(mumumu) - 91.2| < 10
        ## 6. No b-jet
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
        
            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.): return
                if not bjets.size() == 0: return
                self.fillCutflow(CutStage.Final, "WZ3Mu", weight, "Central")
                return "WZ3Mu"
            else:
                if not abs(pair1.M() - 91.2) > 10.: return
                if not abs(pair2.M() - 91.2) > 10.: return
                if not abs((mu1+mu2+mu3).M() - 91.2) < 10.: return
                if not bjets.size() == 0: return
                self.fillCutflow(CutStage.Final, "ZG3Mu", weight, "Central")
                return "ZG3Mu"
        return

    def getConvLepton(self, muons: RVec[Muon], electrons: RVec[Electron]) -> Lepton:
        if electrons.size() == 1 and muons.size() == 2:
            return electrons.at(0);
        elif muons.size() == 3:
            mu1, mu2, mu3 = tuple(muons)
            if mu1.Charge() == mu2.Charge():
                pair1, pair2 = mu1+mu3, mu2+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu2
            elif mu1.Charge() == mu3.Charge():
                pair1, pair2 = mu1+mu2, mu2+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1, pair2 = mu1+mu2, mu1+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu2
                else:
                    return mu3
        else:
            raise NotImplementedError(f"Wrong number of muons (muons.size())")

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
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), totWeight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), totWeight, 100, 0., 100.)

        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), totWeight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), totWeight, 100, 0., 100.)
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), totWeight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), totWeight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), totWeight, 15, 0., 15.)
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), totWeight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), totWeight, 64, -3.2, 3.2)

        ## fill ZCand and conversion candidate
        if "ZG" in channel:
            ZCand = None
            if "1E2Mu" in channel: ZCand = electrons.at(0) + muons.at(0) + muons.at(1)
            if "3Mu" in channel:   ZCand = muons.at(0) + muons.at(1) + muons.at(2)
            convLep = self.getConvLepton(muons, electrons)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/convLep/pt", convLep.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/convLep/eta", convLep.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/convLep/phi", convLep.Phi(), totWeight, 64, -3.2, 3.2)
        if "WZ" in channel:
            ZCand = None
            if "1E2Mu" in channel: 
                ZCand = muons.at(0) + muons.at(1)
                nZCand = muons.at(0) + muons.at(1)
                lep = electrons.at(0)
            if "3Mu" in channel:
                mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
                pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
                if abs(pair1.M() - 91.2) < abs(pair2.M() - 91.2):
                    ZCand, nZCand = pair1, pair2
                    lep = mu_ss2
                else:
                    ZCand, nZCand = pair2, pair1
                    lep = mu_ss1
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), totWeight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), totWeight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), totWeight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/lep/pt", lep.Pt(), totWeight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/lep/eta", lep.Eta(), totWeight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/lep/phi", lep.Phi(), totWeight, 64, -3.2, 3.2)
