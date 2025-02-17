#ifndef GeneratorStudy_h
#define GeneratorStudy_h

#include "AnalyzerCore.h"

class GeneratorStudy : public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();

    RVec<LHE> LHEs;
    RVec<GenJet> GenJets;
    RVec<GenDressedLepton> GenDressedLeptons;

    GeneratorStudy();
    ~GeneratorStudy();
};

#endif
