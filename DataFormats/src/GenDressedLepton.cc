#include "GenDressedLepton.h"

ClassImp(GenDressedLepton)

GenDressedLepton::GenDressedLepton() : Particle() {
    j_pdgId = 0;
    j_hasTauAnc = false;
}

GenDressedLepton::~GenDressedLepton() {}

//////////