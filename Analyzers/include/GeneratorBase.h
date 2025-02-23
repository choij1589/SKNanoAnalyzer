#ifndef GeneratorBase_h
#define GeneratorBase_h

#include "AnalyzerCore.h"

class GeneratorBase : public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();

    RVec<LHE> LHEObjects;
    RVec<Gen> GenObjects;
    RVec<GenJet> GenJetObjects;
    RVec<GenDressedLepton> GenDressedLeptonObjects;
    RVec<GenIsolatedPhoton> GenIsolatedPhotonObjects;
    RVec<GenVisTau> GenVisTauObjects;

    GeneratorBase();
    ~GeneratorBase();
};

#endif