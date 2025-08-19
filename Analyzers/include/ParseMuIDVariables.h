#ifndef ParseMuIDVariables_h
#define ParseMuIDVariables_h

#include "AnalyzerCore.h"

class ParseMuIDVariables: public AnalyzerCore {
private:
    TTree *Events;
    // Event variables
    float genWeight;
    unsigned int nMuons;
    float pt[20];
    float eta[20];
    int lepType[20];
    int nearestJetFlavour[20];
    // ID variables
    float isPOGMediumId[20];
    float dZ[20];
    float sip3d[20];
    float tkRelIso[20];
    float miniPFRelIso[20];

    // Trigger matching variables
    bool isEMuTrigMatched[20];
    bool isIsoMuTrigMatched[20];
    
    // For event selection
    //RVec<TString> EMuTriggers;
    RVec<TString> SglElTriggers;

public:
    ParseMuIDVariables();
    ~ParseMuIDVariables();

    void initializeAnalyzer();
    void executeEvent();
    void WriteHist();
    bool PassSLT(const Electron &el, const RVec<TrigObj> &trigObjs);
    bool PassEMT(const Muon &mu, const RVec<TrigObj> &trigObjs);
    bool PassIsoMuT(const Muon &mu, const RVec<TrigObj> &trigObjs);
};

#endif
