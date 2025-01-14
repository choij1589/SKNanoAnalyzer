#include <algorithm>
// Note that this is a template for MC, for data, it will automatically change uncessary lines in SKNano.py

void [jobname]() {
    [analyzer] module;
    module.SetTreeName("Events");
    module.LogEvery = 1000;
    module.IsDATA = false;
    module.MCSample = "[sample]";
    module.xsec = [xsec];
    module.sumW = [sumW];
    module.sumSign = [sumSign];
    module.SetEra("[era]");
[USERFLAGS]
[SAMPLEPATHS]
[MAXEVENT]
    module.SetOutfilePath("[output]");
    module.Init();
    module.initializeAnalyzer();
    module.Loop();
    module.WriteHist();
}