#ifndef Skim_TriLep_h
#define Skim_TriLep_h

#include "AnalyzerCore.h"

class Skim_TriLep: public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();

    Skim_TriLep();
    ~Skim_TriLep();

    TTree *newtree;

    RVec<TString> triggers;
    void WriteHist();

};

#endif
