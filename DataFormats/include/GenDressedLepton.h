#ifndef GenDressedLepton_h
#define GenDressedLepton_h

#include "Particle.h"

class GenDressedLepton: public Particle {
public:
    GenDressedLepton();
    ~GenDressedLepton();

    inline void SetPdgId(Int_t pdgId) { j_pdgId = pdgId; }
    inline void SetHasTauAnc(Bool_t hasTauAnc) { j_hasTauAnc = hasTauAnc; }

    inline Int_t PdgId() const { return j_pdgId; }
    inline Bool_t HasTauAnc() const { return j_hasTauAnc; }

private:
    Int_t j_pdgId;
    Bool_t j_hasTauAnc;
    ClassDef(GenDressedLepton, 1)
};

#endif