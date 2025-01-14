#include "Gen.h"

ClassImp(Gen)

Gen::Gen() : Particle() {
    j_IsEmpty = true;
    j_index = -1;
    j_PID = 0;
    j_status = 0;
    j_mother_index = -1;
    j_status_flags = -1;   
}

Gen::~Gen() {}

void Gen::SetIsEmpty(bool b){
    j_IsEmpty = b;
}

void Gen::SetIndexPIDStatus(int i, int p, int t){
    j_index = i;
    j_PID = p;
    j_status = t;
}

void Gen::SetMother(int mindex){
    j_mother_index = mindex;
}
void Gen::SetGenStatusFlags(unsigned short f){
    j_status_flags = f;
}

void Gen::Print() const{
    cout << "(Index, PID, Status, MotherIndex, Pt, Eta, Phi, M) = " << Index() << "\t" << PID() << "\t" << Status() << "\t" << MotherIndex() << "\t" << Pt() << "\t" << Eta() << "\t" << Phi() << "\t" << M() << endl;
}
