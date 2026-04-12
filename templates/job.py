import os, sys
sys.path.append(f"{os.getenv("SKNANO_HOME")}/PyAnalyzers")
from ROOT import TString
from ROOT.VecOps import RVec
from [Analyzer] import [Analyzer]

if __name__ == "__main__":
    module = [Analyzer]()
    module.SetTreeName("Events")
    module.LogEvery = 5000
    module.IsDATA = False
    module.MCSample = "[sample]"
    module.xsec = [xsec]
    module.sumW = [sumW]
    module.sumSign = [sumSign]
    module.filterEff = [filterEff]
    module.SetEra("[era]")
    module.SetPeriod("[period]")
[USERFLAGS]
[SAMPLEPATHS]
[MAXEVENT]
    module.SetOutfilePath("[output]")
    module.Init()
    module.initializePyAnalyzer()
    module.Loop()
    module.WriteHist()
