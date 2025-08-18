#ifndef TrigObj_h
#define TrigObj_h

#include "Particle.h"

class TrigObj : public Particle{
public:
    TrigObj();
    ~TrigObj();
    
    // Setting functions
    inline void SetId(int id) { j_id = id; }
    inline void SetFilterBits(int filterBits) { j_filterBits = filterBits; }
    inline void SetRun(int run) { j_run = run; }
    
    // Getter functions
    inline int id() const { return j_id; }
    inline int filterBits() const { return j_filterBits; }
    inline int run() const { return j_run; }
    
    // Helper functions to check object type based on ID
    inline bool isElectron() const { return j_id == 11; }
    inline bool isPhoton() const { return j_id == 22; }
    inline bool isMuon() const { return j_id == 13; }
    inline bool isTau() const { return j_id == 15; }
    inline bool isJet() const { return j_id == 1; }
    inline bool isFatJet() const { return j_id == 6; }
    inline bool isMET() const { return j_id == 2; }
    inline bool isHT() const { return j_id == 3; }
    inline bool isMHT() const { return j_id == 4; }
    inline bool isBoostedTau() const { return j_id == 1515; }
    
    // Helper function to check if a specific bit is set
    inline bool hasBit(int bit) const { return (j_filterBits & (1 << bit)) != 0; }
    
private:
    // Object ID: 11=Electron, 22=Photon, 13=Muon, 15=Tau, 1=Jet, 6=FatJet, 2=MET, 3=HT, 4=MHT, 1515=BoostedTau
    int j_id;
    // Filter bits with associated trigger information
    int j_filterBits;
    // Run number (2 or 3) to determine which filterBits definition to use
    int j_run;
    
    ClassDef(TrigObj, 1)
};

#endif
