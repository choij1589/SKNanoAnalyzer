#ifndef ParseEleIDVariables_h
#define ParseEleIDVariables_h

#include "AnalyzerCore.h"

class ParseEleIDVariables: public AnalyzerCore {
private:
    TTree *Events;
    // Event variables
    float genWeight;
    unsigned int nElectrons;
    float pt[20];
    float scEta[20];
    int lepType[20];
    // Electrons
    // Trigger variables
    float sieie[20];
    float deltaEtaInSC[20];
    float deltaPhiInSeed[20];
    float hoe[20];
    float ecalPFClusterIso[20];
    float hcalPFClusterIso[20];
    float rho[20];
    float dr03TkSumPt[20];
    
    // ID variables
    bool isMVANoIsoWP90[20];
    bool isPOGMedium[20];
    bool isPOGTight[20];
    bool convVeto[20];
    int  lostHits[20];
    float dZ[20];
    float sip3d[20];
    float miniPFRelIso[20];
    float mvaNoIso[20];
    int   nearestJetFlavour[20];

    // Trigger matching variables
    bool isIsoElTrigMatched[20];
    bool isEMuTrigMatched[20];

    // For event selection
    RVec<TString> SglMuTriggers;
    //RVec<TString> EMuTriggers;

public:
    ParseEleIDVariables();
    ~ParseEleIDVariables();

    void initializeAnalyzer();
    void executeEvent();
    void WriteHist();
    bool PassSLT(const Muon &mu, const RVec<TrigObj> &trigObjs);
    bool PassIsoElT(const Electron &el, const RVec<TrigObj> &trigObjs);
    bool PassEMT(const Electron &el, const RVec<TrigObj> &trigObjs);
};

#endif
