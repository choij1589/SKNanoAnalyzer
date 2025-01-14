#ifndef GenAnalyzer_h
#define GenAnalyzer_h

#include "AnalyzerCore.h"

class GenAnalyzer : public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();

    RVec<Gen> AllGens;
    RVec<LHE> AllLHEs;

    GenAnalyzer();
    ~GenAnalyzer();
};

#endif
