#ifndef ClosDiLepTrigs_h
#define ClosDiLepTrigs_h

#include "DiLeptonBase.h"

class ClosDiLepTrigs: public DiLeptonBase {
public:
    ClosDiLepTrigs();
    ~ClosDiLepTrigs();

    enum class Channel {
        NONE,
        DIMU,
        EMU,
        EMUMU,
        MUMUMU
    };

    struct RecoObjects {
        RVec<Muon> looseMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> looseElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> tightJets;
        RVec<Jet> tightJets_vetoLep;
        RVec<Jet> bjets;
        Particle METv;
    };

    void initializeAnalyzer();
    void executeEvent();

private:
    Channel channel;

    RecoObjects defineObjects(Event& ev, 
                             const RVec<Muon>& rawMuons, 
                             const RVec<Electron>& rawElectrons, 
                             const RVec<Jet>& rawJets);
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects);
    float getWeight(const Channel& channel, const Event& event, 
                   const RecoObjects& recoObjects);
    void fillObjects(const Channel& channel, const RecoObjects& recoObjects, 
                    float weight);

    TString channelToString(const Channel& channel);
};

#endif
