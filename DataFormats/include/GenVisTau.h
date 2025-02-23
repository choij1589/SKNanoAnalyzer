#ifndef GenVisTau_h
#define GenVisTau_h

#include "Particle.h"

class GenVisTau: public Particle {
public:
    GenVisTau();
    ~GenVisTau();

    inline void SetCharge(Int_t charge) { j_charge = charge; }
    inline void SetGenPartIdxMother(Int_t genPartIdxMother) { j_genPartIdxMother = genPartIdxMother; }
    inline void SetStatus(Int_t status) { j_status = status; }

    inline Int_t Charge() const { return j_charge; }
    inline Int_t GenPartIdxMother() const { return j_genPartIdxMother; }
    inline Int_t Status() const { return j_status; }

private:
    Int_t j_charge;
    Int_t j_genPartIdxMother;
    Int_t j_status;

    ClassDef(GenVisTau, 1)
};


#endif