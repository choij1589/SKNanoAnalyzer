#include "Electron.h"

ClassImp(Electron)

Electron::Electron() {
    this->SetLeptonFlavour(ELECTRON);
    // uncertainties
    j_energyErr = -999.;
    // ID variables
    j_convVeto = false;
    j_lostHits = 0;
    j_seedGain = 0;
    j_tightCharge = 0;
    j_sieie = -999.;
    j_hoe = -999.;
    j_eInvMinusPInv = -999.;
    j_dr03EcalRecHitSumEt = -999.;
    j_dr03HcalDepth1TowerSumEt = -999.;
    j_dr03TkSumPt = -999.;
    j_dr03TkSumPtHEEP = -999.;
    // supercluster
    j_deltaEtaSC = -999.;
    // ID bits
    j_mvaIso_WP80 = false;
    j_mvaIso_WP90 = false;
    j_mvaNoIso_WP80 = false;
    j_mvaNoIso_WP90 = false;
    j_cutBased_HEEP = false;
    j_cutBased = 0;
    j_mvaIso = -999.;
    j_mvaNoIso = -999.;
    j_mvaTTH = -999.;
    // others
    j_r9 = -999.;
    j_genPartFlav = 0;
    j_genPartIdx = -1;
}

Electron::~Electron() {}

void Electron::SetBIDBit(BooleanID id, bool idbit) {
    switch (id) {
        case BooleanID::MVAISOWP80:            j_mvaIso_WP80 = idbit; break;
        case BooleanID::MVAISOWP90:            j_mvaIso_WP90 = idbit; break;
        case BooleanID::MVANOISOWP80:          j_mvaNoIso_WP80 = idbit; break;
        case BooleanID::MVANOISOWP90:          j_mvaNoIso_WP90 = idbit; break;
        case BooleanID::CUTBASEDHEEP:          j_cutBased_HEEP = idbit; break;
        default: break;
    }
}

void Electron::SetCBIDBit(CutBasedID id, unsigned int idbit) {
    switch (id) {
        case CutBasedID::CUTBASED:             j_cutBased = idbit; break;
        default: break;
    }
}

void Electron::SetMVA(MVATYPE type, float score) {
    switch (type) {
        case MVATYPE::MVAISO:                    j_mvaIso = score; break;
        case MVATYPE::MVANOISO:                  j_mvaNoIso = score; break;
        case MVATYPE::MVATTH:                    j_mvaTTH = score; break;
        default: break;
    }
}

bool Electron::PassID(const TString ID) const {
    // always veto gap
    // cout << "Eta: " << Eta() << "\t" << scEta() << endl;
    if (etaRegion() == ETAREGION::GAP) return false;

    // POG
    if (ID == "POGVeto")          return (int)(CutBased()) >= (int)(WORKINGPOINT::VETO);
    if (ID == "POGLoose")         return (int)(CutBased()) >= (int)(WORKINGPOINT::LOOSE);
    if (ID == "POGMedium")        return (int)(CutBased()) >= (int)(WORKINGPOINT::MEDIUM);
    if (ID == "POGTight")         return (int)(CutBased()) >= (int)(WORKINGPOINT::TIGHT);
    if (ID == "POGHEEP")          return isCutBasedHEEP();
    if (ID == "POGMVAIsoWP80")    return isMVAIsoWP80();
    if (ID == "POGMVAIsoWP90")    return isMVAIsoWP90();
    if (ID == "POGMVANoIsoWP80")  return isMVANoIsoWP80();
    if (ID == "POGMVANoIsoWP90")  return isMVANoIsoWP90();
    // if (ID == "TopHNTight_16a")   return Pass_TopHNTight_16a();
    // if (ID == "TopHNTight_16b")   return Pass_TopHNTight_16b();
    // if (ID == "TopHNTight_17")    return Pass_TopHNTight_17();
    // if (ID == "TopHNTight_18")    return Pass_TopHNTight_18();
    // if (ID == "TopHNLoose_16a")   return Pass_TopHNLoose_16a();
    // if (ID == "TopHNLoose_16b")   return Pass_TopHNLoose_16b();
    // if (ID == "TopHNLoose_17")    return Pass_TopHNLoose_17();
    //if (ID == "TopHNLoose_18")    return Pass_TopHNLoose_18();

    cerr << "[Electron::PassID] " << ID << " is not implemented" << endl;
    exit(ENODATA);

    return false;
}

bool Electron::PassID(ElectronID ID) const {
    switch (ID) {
        case ElectronID::POG_VETO:         return PassID("POGVeto");
        case ElectronID::POG_LOOSE:        return PassID("POGLoose");
        case ElectronID::POG_MEDIUM:       return PassID("POGMedium");
        case ElectronID::POG_TIGHT:        return PassID("POGTight");
        case ElectronID::POG_HEEP:         return PassID("POGHEEP");
        case ElectronID::POG_MVAISO_WP80:      return PassID("POGMVAIsoWP80");
        case ElectronID::POG_MVAISO_WP90:      return PassID("POGMVAIsoWP90");
        case ElectronID::POG_MVANOISO_WP80:    return PassID("POGMVANoIsoWP80");
        case ElectronID::POG_MVANOISO_WP90:    return PassID("POGMVANoIsoWP90");
        default: break;
    }
    return false;
}