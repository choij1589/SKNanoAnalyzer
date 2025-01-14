#include "LHE.h"

ClassImp(LHE)

LHE::LHE() : Particle()
{
    j_incomingPz = -999.;
    j_spin = -999;
    j_status = -999;
}

void LHE::Print() const {
    cout << "LHE: " << Pt() << ", " << Eta() << ", " << Phi() << ", " << M() << endl;
}

LHE::~LHE() {}
