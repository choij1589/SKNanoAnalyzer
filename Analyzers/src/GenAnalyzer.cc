#include "GenAnalyzer.h"

GenAnalyzer::GenAnalyzer() {}
GenAnalyzer::~GenAnalyzer() {}

void GenAnalyzer::initializeAnalyzer() {
}
void GenAnalyzer::executeEvent() {
    AllGens = GetAllGens();
    AllLHEs = GetAllLHEs();

    //cout << "AllGens: " << AllGens.size() << endl;
    //cout << "AllLHEs: " << AllLHEs.size() << endl;

    for (const auto &gen: AllGens) {
        gen.Print();
    }
    for (const auto &lhe: AllLHEs) {
        lhe.Print();
    }
}
