import os
from ROOT import TString
from ROOT.VecOps import RVec
from MatrixSelector import MatrixSelector
from ROOT import JetTagging
from ROOT import Event, Lepton, Muon, Electron, Jet

class CRMatrixSelector(MatrixSelector):
    def __init__(self):
        super().__init__()

    def selectEvent(self, ev: Event,
                          recoObjects: dict) -> str:
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
        
        #### Not all leptons tight
        if self.Run1E2Mu:
            if not is1E2Mu: return
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return

        if self.Run3Mu:
            if not is3Mu: return
            if tightMuons.size() == looseMuons.size(): return

        # 1E2Mu ZGamma
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. |M(mumue) - 91.2| < 10
        ## 5. No b-jet
        if self.Run1E2Mu:
            if not ev.PassTrigger(self.EMuTriggers): return
                
            mu1, mu2 = looseMuons.at(0), looseMuons.at(1)
            ele = looseElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return
            if not mu1.Charge()+mu2.Charge() == 0: return
            pair = mu1 + mu2
            if not (pair.M() > 12.): return
            
            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair.M() - 91.2) < 10.): return
                if not bjets.size() == 0: return
                return "WZ1E2Mu"
            else:
                if not (abs(pair.M() - 91.2) > 10.): return
                if not abs((mu1+mu2+ele).M() - 91.2) < 10.: return
                if not bjets.size() == 0: return
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
                
            mu1, mu2, mu3 = tuple(looseMuons)
            if not mu1.Pt() > 20.: return
            if not mu2.Pt() > 10.: return
            if not mu3.Pt() > 10.: return
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(looseMuons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            if not pair1.M() > 12.: return
            if not pair2.M() > 12.: return
        
            if METv.Pt() > 40.:
                # WZ control region
                if not (abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.): return
                if not bjets.size() == 0: return
                return "WZ3Mu"
            else:
                if not abs(pair1.M() - 91.2) > 10.: return
                if not abs(pair2.M() - 91.2) > 10.: return
                if not abs((mu1+mu2+mu3).M() - 91.2) < 10.: return
                if not bjets.size() == 0: return
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
                          weight: float,
                          syst: str = "Central"):
        muons = recoObjects["looseMuons"]
        electrons = recoObjects["looseElectrons"]
        jets = recoObjects["jets"]
        bjets = recoObjects["bjets"]
        METv = recoObjects["METv"]
        
        ## fill base observables
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)

        # fill ZCand and conversion candidate
        if "ZG" in channel:
            ZCand = None
            if "1E2Mu" in channel:
                ZCand = electrons.at(0) + muons.at(0) + muons.at(1)
            if "3Mu" in channel:
                ZCand = muons.at(0) + muons.at(1) + muons.at(2)
            convLep = self.getConvLepton(muons, electrons)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/convLep/pt", convLep.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/convLep/eta", convLep.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/convLep/phi", convLep.Phi(), weight, 64, -3.2, 3.2)
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
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/lep/pt", lep.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/lep/eta", lep.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/lep/phi", lep.Phi(), weight, 64, -3.2, 3.2)
