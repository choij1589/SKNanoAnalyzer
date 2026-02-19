#ifndef GetEffLumi_h 
#define GetEffLumi_h

#include "AnalyzerCore.h"

class GetEffLumi : public AnalyzerCore {
public:
    void executeEvent();
    void WriteHist();
    GetEffLumi();
    ~GetEffLumi();
    void initializeAnalyzer();

private:
    long long genFilterPassed = 0;
    long long genFilterTotal = 0;
};

#endif

