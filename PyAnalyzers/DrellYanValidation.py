## Note for writing python analyzers
## 1. Inherit from Base CPP Analyzer - It should have initializeAnalyzer() but executeEvent() is not mandatory
## 2. Do not use the same name for the analyzer class in python and cpp - it will cause recursive error
from ROOT import gSystem
from ROOT import GeneratorBase
from ROOT import TString, TMath
from ROOT.VecOps import RVec
from ROOT import LHE, Gen, GenJet, GenDressedLepton, GenIsolatedPhoton, GenVisTau

class GenValidation(GeneratorBase):
    def __init__(self):
        super().__init__()

    def initializeFromBase(self):
        self.initializeAnalyzer()
        self.genObjects = {}
    
    def executeEvent(self):
        genRawObjects = self.reNewGenRawObjects()
        self.selectLHEObjects(genRawObjects)
        self.selectGenObjects(genRawObjects)
        self.selectGenJetObjects(genRawObjects)
        #print(self.genObjects["GEN-leptons"])

        self.fillObjects()

    def reNewGenRawObjects(self):
        LHEObjects = self.GetAllLHEs()
        GenObjects = self.GetAllGens()
        GenJetObjects = self.GetAllGenJets()
        GenDressedLeptonObjects = self.GetAllGenDressedLeptons()
        GenIsolatedPhotonObjects = self.GetAllGenIsolatedPhotons()
        GenVisTauObjects = self.GetAllGenVisTaus()

        genRawObjects = {
            "LHE": LHEObjects,
            "Gen": GenObjects,
            "GenJet": GenJetObjects,
            "GenDressedLepton": GenDressedLeptonObjects,
            "GenIsolatedPhoton": GenIsolatedPhotonObjects,
            "GenVisTau": GenVisTauObjects
        }
        return genRawObjects

    def selectLHEObjects(self, genRawObjects):
        LHEObjects = genRawObjects["LHE"]
        lhe_leptons = RVec("LHE")()
        lhe_outgoing_jets = RVec("LHE")()

        for lhe in LHEObjects:
            if lhe.Status() == -1: continue # incoming particles
            pdgId = abs(lhe.PdgId())
            if pdgId in [11, 13, 15]:
                lhe_leptons.emplace_back(lhe)
            else:
                lhe_outgoing_jets.emplace_back(lhe)
        
        # Define channel based on LHE leptons and outgoing jets
        assert lhe_leptons.size() == 2, f"Unexpected number of leptons: {lhe_leptons.size()}"
        channel = f"DY{lhe_outgoing_jets.size()}"
        mapping = {(11, 11): "_ee", (13, 13): "_mumu", (15, 15): "_tautau"}
        key = tuple(sorted([abs(lhe_leptons[0].PdgId()), abs(lhe_leptons[1].PdgId())]))
        try:
            channel += mapping[key]
        except KeyError:
            raise ValueError(f"Unexpected leptons: {lhe_leptons[0].PdgId()}, {lhe_leptons[1].PdgId()}")        
        
        self.genObjects["channel"] = channel
        self.genObjects["LHE-leptons"] = lhe_leptons
        self.genObjects["LHE-outgoing_jets"] = lhe_outgoing_jets

    def selectGenObjects(self, genRawObjects):
        # Find Gen leptons that can be matched to LHE leptons
        LHE_leptons = self.genObjects["LHE-leptons"]
        GenObjects = genRawObjects["Gen"]

        gen_leptons = RVec("Gen")()
        # Using GetLeptonType()
        # Checked event level matching eff. > 99.9% for electrons and muions
        # For tau leptons, most of the events are missing due to hadronically decaying taus
        for gen in GenObjects:
            lepton_type = self.GetLeptonType(gen, GenObjects)
            if lepton_type in [1, 3]:    # Is EWPrompt or EWtau daughter
                gen_leptons.emplace_back(gen)
        self.genObjects["GEN-leptons"] = gen_leptons

    def selectGenJetObjects(self, genRawObjects):
        GenJetObjects = genRawObjects["GenJet"]
        gen_leptons = self.genObjects["GEN-leptons"]
        gen_jets = RVec("GenJet")()

        for gen_jet in GenJetObjects:
            # Clean jets overlapping with leptons
            if any(gen_jet.DeltaR(gen_lepton) < 0.4 for gen_lepton in gen_leptons): continue
            gen_jets.emplace_back(gen_jet)
        self.genObjects["GEN-jets"] = gen_jets

    def fillObjects(self):
        channel = self.genObjects["channel"]
        LHE_leptons = self.genObjects["LHE-leptons"]
        LHE_jets = self.genObjects["LHE-outgoing_jets"]
        GEN_leptons = self.genObjects["GEN-leptons"]
        GEN_jets = self.genObjects["GEN-jets"]

        # Fill LHE distributions
        LHE_ZCand = LHE_leptons[0] + LHE_leptons[1]
        for idx, jet in enumerate(LHE_jets, start=1):
            self.FillHist(f"Inclusive/LHE/jets/{idx}/pt", jet.Pt(), 1, 100, 0, 100)
            self.FillHist(f"Inclusive/LHE/jets/{idx}/eta", jet.Eta(), 1, 200, -10, 10)
            self.FillHist(f"Inclusive/LHE/jets/{idx}/phi", jet.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"Inclusive/LHE/jets/{idx}/mass", jet.M(), 1, 100, 0, 100)
        self.FillHist("Inclusive/LHE/jets/size", LHE_jets.size(), 1, 10, 0, 10)

        self.FillHist("Inclusive/LHE/Z_pt", LHE_ZCand.Pt(), 1, 100, 0, 100)
        self.FillHist("Inclusive/LHE/Z_eta", LHE_ZCand.Eta(), 1, 200, -10, 10)
        self.FillHist("Inclusive/LHE/Z_phi", LHE_ZCand.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
        self.FillHist("Inclusive/LHE/Z_mass", LHE_ZCand.M(), 1, 300, 0, 300)

        for idx, lepton in enumerate(LHE_leptons, start=1):
            self.FillHist(f"Inclusive/LHE/leptons/{idx}/pt", lepton.Pt(), 1, 100, 0, 100)
            self.FillHist(f"Inclusive/LHE/leptons/{idx}/eta", lepton.Eta(), 1, 200, -10, 10)
            self.FillHist(f"Inclusive/LHE/leptons/{idx}/phi", lepton.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"Inclusive/LHE/leptons/{idx}/mass", lepton.M(), 1, 100, 0, 100)

        for idx, jet in enumerate(LHE_jets, start=1):
            self.FillHist(f"{channel}/LHE/jets/{idx}/pt", jet.Pt(), 1, 100, 0, 100)
            self.FillHist(f"{channel}/LHE/jets/{idx}/eta", jet.Eta(), 1, 200, -10, 10)
            self.FillHist(f"{channel}/LHE/jets/{idx}/phi", jet.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"{channel}/LHE/jets/{idx}/mass", jet.M(), 1, 100, 0, 100)

        self.FillHist(f"{channel}/LHE/Z_pt", LHE_ZCand.Pt(), 1, 100, 0, 100)
        self.FillHist(f"{channel}/LHE/Z_eta", LHE_ZCand.Eta(), 1, 200, -10, 10)
        self.FillHist(f"{channel}/LHE/Z_phi", LHE_ZCand.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
        self.FillHist(f"{channel}/LHE/Z_mass", LHE_ZCand.M(), 1, 300, 0, 300)
        
        for idx, lepton in enumerate(LHE_leptons, start=1):
            self.FillHist(f"{channel}/LHE/leptons/{idx}/pt", lepton.Pt(), 1, 100, 0, 100)
            self.FillHist(f"{channel}/LHE/leptons/{idx}/eta", lepton.Eta(), 1, 200, -10, 10)
            self.FillHist(f"{channel}/LHE/leptons/{idx}/phi", lepton.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"{channel}/LHE/leptons/{idx}/mass", lepton.M(), 1, 100, 0, 100)

        # Fill GEN distributions
        if not (GEN_leptons.size() == 2): return
        GEN_ZCand = GEN_leptons[0] + GEN_leptons[1]
        for idx, jet in enumerate(GEN_jets, start=1):
            self.FillHist(f"Inclusive/GEN/jets/{idx}/pt", jet.Pt(), 1, 100, 0, 100)
            self.FillHist(f"Inclusive/GEN/jets/{idx}/eta", jet.Eta(), 1, 200, -10, 10)
            self.FillHist(f"Inclusive/GEN/jets/{idx}/phi", jet.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"Inclusive/GEN/jets/{idx}/mass", jet.M(), 1, 100, 0, 100)
        self.FillHist("Inclusive/GEN/jets/size", GEN_jets.size(), 1, 10, 0, 10)
        
        self.FillHist(f"Inclusive/GEN/Z_pt", GEN_ZCand.Pt(), 1, 100, 0, 100)
        self.FillHist(f"Inclusive/GEN/Z_eta", GEN_ZCand.Eta(), 1, 200, -10, 10)
        self.FillHist(f"Inclusive/GEN/Z_phi", GEN_ZCand.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
        self.FillHist(f"Inclusive/GEN/Z_mass", GEN_ZCand.M(), 1, 300, 0, 300)

        for idx, lepton in enumerate(GEN_leptons, start=1):
            self.FillHist(f"Inclusive/GEN/leptons/{idx}/pt", lepton.Pt(), 1, 100, 0, 100)
            self.FillHist(f"Inclusive/GEN/leptons/{idx}/eta", lepton.Eta(), 1, 200, -10, 10)
            self.FillHist(f"Inclusive/GEN/leptons/{idx}/phi", lepton.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"Inclusive/GEN/leptons/{idx}/mass", lepton.M(), 1, 100, 0, 100)

        for idx, jet in enumerate(GEN_jets, start=1):
            self.FillHist(f"{channel}/GEN/jets/{idx}/pt", jet.Pt(), 1, 100, 0, 100)
            self.FillHist(f"{channel}/GEN/jets/{idx}/eta", jet.Eta(), 1, 200, -10, 10)
            self.FillHist(f"{channel}/GEN/jets/{idx}/phi", jet.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"{channel}/GEN/jets/{idx}/mass", jet.M(), 1, 100, 0, 100)
        
        self.FillHist(f"{channel}/GEN/Z_pt", GEN_ZCand.Pt(), 1, 100, 0, 100)
        self.FillHist(f"{channel}/GEN/Z_eta", GEN_ZCand.Eta(), 1, 200, -10, 10)
        self.FillHist(f"{channel}/GEN/Z_phi", GEN_ZCand.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
        self.FillHist(f"{channel}/GEN/Z_mass", GEN_ZCand.M(), 1, 300, 0, 300)
        
        for idx, lepton in enumerate(GEN_leptons, start=1):
            self.FillHist(f"{channel}/GEN/leptons/{idx}/pt", lepton.Pt(), 1, 100, 0, 100)
            self.FillHist(f"{channel}/GEN/leptons/{idx}/eta", lepton.Eta(), 1, 200, -10, 10)
            self.FillHist(f"{channel}/GEN/leptons/{idx}/phi", lepton.Phi(), 1, 72, -TMath.Pi(), TMath.Pi())
            self.FillHist(f"{channel}/GEN/leptons/{idx}/mass", lepton.M(), 1, 100, 0, 100)

if __name__ == "__main__":
    module = GenValidation()
    module.SetTreeName("Events")
    module.LogEvery = 1000
    module.IsDATA = False
    module.MCSample = "DYJets"
    module.SetEra("2023")
    if not module.AddFile("/Users/choij/Sync/workspace/SKNanoAnalyzer/test/NANOGEN/MG4GPU/DY012j_CPPAVX2/NANOGEN.root"): exit(1)
    module.MaxEvent = int(module.fChain.GetEntries())
    module.SetOutfilePath("test.root")
    module.Init()
    module.initializeFromBase()
    module.Loop()
    module.WriteHist()
