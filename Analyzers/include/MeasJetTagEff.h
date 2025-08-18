#ifndef MeasJetTagEff_h
#define MeasJetTagEff_h

#include "DiLeptonBase.h"

class MeasJetTagEff : public DiLeptonBase {
public:
    MeasJetTagEff();
    ~MeasJetTagEff();
    void initializeAnalyzer();
    void executeEvent();

    RVec<Muon> looseMuons, tightMuons;
    RVec<Electron> looseElectrons, tightElectrons;
    RVec<Jet> tightJets, tightJets_vetoLep;

    vector<float> vec_etabins;
    vector<float> vec_ptbins;

    float PtMax;

    int NEtaBin;
    int NPtBin;

    float *etabins;
    float *ptbins;

    RVec<JetTagging::JetFlavTagger> Taggers;
    RVec<JetTagging::JetFlavTaggerWP> WPs;
};

#endif
