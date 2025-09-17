#ifndef TriggerStudy_h
#define TriggerStudy_h

#include "AnalyzerCore.h"

class TriggerStudy: public AnalyzerCore {
public:
    TriggerStudy();
    ~TriggerStudy();

    // Trigger Setup
    bool RunEMuTrigs;
    bool RunEMuTrigsWithSglMuTrigs;
    bool RunEMuTrigsWithSglElTrigs;
    bool RunEMuTrigsWithDblMuTrigs;

    // IDs
    IDContainer *MuonIDs, *ElectronIDs;
    
    // Triggers
    RVec<TString> EMuTriggers;

    void initializeAnalyzer();
    void executeEvent();

    enum class Channel {
        NONE,
        SR1E2MU,
        SR3MU
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::SR1E2MU) return "SR1E2Mu";
        if (ch == Channel::SR3MU) return "SR3Mu";
        return "None";
    }

    enum class CutStage {
        Initial = 0,
        NoiseFilter = 1,
        LeptonSelection = 2,
        Trigger = 3,
        KinematicCuts = 4,
        OSMuonPair = 5,
        JetRequirements = 6,
        BjetRequirements = 7,
        Final = 8
    };

    struct RecoObjects {
        RVec<Muon> vetoMuons;
        RVec<Muon> tightMuons;
        RVec<Electron> vetoElectrons;
        RVec<Electron> tightElectrons;
        RVec<Jet> tightJets;
        RVec<Jet> bjets;
        Particle METv;
    };
    struct GenObjects {
        RVec<Gen> genParts;
        RVec<GenJet> genJets;
    };

private:
    RecoObjects defineObjects(Event& ev, RVec<Muon>& allMuons, RVec<Electron>& allElectrons, RVec<Jet>& allJets, const TString& syst="Central");
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects, const GenObjects& genObjects, const TString& syst="Central");
    void fillCutflow(CutStage stage, const Channel& channel, float weight, const TString& syst);
    void fillObjects(const Channel& channel, const RecoObjects& recoObjects, const GenObjects& genObjects, float weight, const TString& syst="Central");
};

#endif