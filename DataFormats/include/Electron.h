#ifndef Electron_h
#define Electron_h

#include "TString.h"
#include "Lepton.h"

// Need update
// Missing variables (compared to SKFlat)
// energy / resolution corrections & errors
// ID variables: j_e2x5OverE5x5, j_e1x5OverE5x5, j_ecalPFClusterIso, j_hcalPFClusterIso, j_dr03HcalTowerSumEt
// SuperCluster: phi, E
// others: j_ea

class Electron : public Lepton {
public:
    Electron();
    ~Electron();

    enum class ETAREGION {IB, OB, GAP, EC};
    inline ETAREGION etaRegion() const {
        if (fabs(scEta()) < 0.8)        return ETAREGION::IB;
        else if (fabs(scEta()) < 1.444) return ETAREGION::OB;
        else if (fabs(scEta()) < 1.566) return ETAREGION::GAP;
        else return ETAREGION::EC;
    }

    enum class ElectronID {
        NOCUT,
        POG_VETO,
        POG_LOOSE,
        POG_MEDIUM,
        POG_TIGHT,
        POG_HEEP,
        POG_MVAISO_WP80,
        POG_MVAISO_WP90,
        POG_MVAISO_WPL,
        POG_MVANOISO_WP80,
        POG_MVANOISO_WP90,
        POG_MVANOISO_WPL,
    };

    void SetConvVeto(bool convVeto) { j_convVeto = convVeto; }
    inline bool ConvVeto() const { return j_convVeto; }

    void SetScEta(float scEta) { j_scEta = scEta; }
    inline float scEta() const { return j_scEta; }

    void SetScPhi(float scPhi) { j_scPhi = scPhi; }
    inline float scPhi() const { return j_scPhi; }

    void SetDeltaEtaInSC(float deltaEtaInSC) { j_deltaEtaInSC = deltaEtaInSC; }
    inline float deltaEtaInSC() const { return j_deltaEtaInSC; }

    void SetDeltaPhiInSC(float deltaPhiInSC) { j_deltaPhiInSC = deltaPhiInSC; }
    inline float deltaPhiInSC() const { return j_deltaPhiInSC; }

    void SetDeltaEtaInSeed(float deltaEtaInSeed) { j_deltaEtaInSeed = deltaEtaInSeed; }
    inline float deltaEtaInSeed() const { return j_deltaEtaInSeed; }

    void SetDeltaPhiInSeed(float deltaPhiInSeed) { j_deltaPhiInSeed = deltaPhiInSeed; }
    inline float deltaPhiInSeed() const { return j_deltaPhiInSeed; }
    
    void SetLostHits(unsigned char lostHits) { j_lostHits = lostHits; }
    inline unsigned char LostHits() const { return j_lostHits; }

    void SetPFClusterIso(float ecalPFClusterIso, float hcalPFClusterIso) {
        j_ecalPFClusterIso = ecalPFClusterIso;
        j_hcalPFClusterIso = hcalPFClusterIso;
    }
    inline float ecalPFClusterIso() const { return j_ecalPFClusterIso; }
    inline float hcalPFClusterIso() const { return j_hcalPFClusterIso; }
    void SetSeedGain(unsigned char seedGain) { j_seedGain = seedGain; }
    inline unsigned char SeedGain() const { return j_seedGain; }

    void SetTightCharge(unsigned char tightCharge) { j_tightCharge = tightCharge; }
    inline unsigned char TightCharge() const { return j_tightCharge; } // 0: None, 1: isGsfScPixChargeConsistent, 2: isGsfCtfScPixChargeConsistent

    void SetSieie(float sieie) { j_sieie = sieie; }
    inline float sieie() const { return j_sieie; }

    void SetHoe(float hoe) { j_hoe = hoe; }
    inline float hoe() const { return j_hoe; }

    void SetEInvMinusPInv(float eInvMinusPInv) { j_eInvMinusPInv = eInvMinusPInv; }
    inline float eInvMinusPInv() const { return j_eInvMinusPInv; }

    void SetDr03EcalRecHitSumEt(float dr03EcalRecHitSumEt) { j_dr03EcalRecHitSumEt = dr03EcalRecHitSumEt; }
    inline float dr03EcalRecHitSumEt() const { return j_dr03EcalRecHitSumEt; }

    void SetDr03HcalDepth1TowerSumEt(float dr03HcalDepth1TowerSumEt) { j_dr03HcalDepth1TowerSumEt = dr03HcalDepth1TowerSumEt; }
    inline float dr03HcalDepth1TowerSumEt() const { return j_dr03HcalDepth1TowerSumEt; }

    void SetDr03TkSumPt(float dr03TkSumPt) { j_dr03TkSumPt = dr03TkSumPt; }
    inline float dr03TkSumPt() const { return j_dr03TkSumPt; }

    void SetDr03TkSumPtHEEP(float dr03TkSumPtHEEP) { j_dr03TkSumPtHEEP = dr03TkSumPtHEEP; }
    inline float dr03TkSumPtHEEP() const { return j_dr03TkSumPtHEEP; }

    void SetR9(float r9) { j_r9 = r9; }
    inline float r9() const { return j_r9; }

    void SetGenPartFlav(unsigned char genPartFlav) { j_genPartFlav = genPartFlav; }
    inline unsigned char GenPartFlav() const { return j_genPartFlav; }

    void SetGenPartIdx(short genPartIdx) { j_genPartIdx = genPartIdx; }
    inline short GenPartIdx() const { return j_genPartIdx; }

    void SetJetIdx(short jetIdx) { j_jetIdx = jetIdx; }
    inline short JetIdx() const { return j_jetIdx; }

    void SetRho(float rho) { j_rho = rho; }
    inline float rho() const { return j_rho; }

    void SetEnergyErr(float energyErr) { j_energyErr = energyErr; }
    inline float energyErr() const { return j_energyErr; }

    void SetEnergyResUnc(float dEsigmaUp, float dEsigmaDown) {
        j_dEsigmaUp = dEsigmaUp;
        j_dEsigmaDown = dEsigmaDown;
    }
    inline float dEsigmaUp() const { return j_dEsigmaUp; }
    inline float dEsigmaDown() const { return j_dEsigmaDown; }

    // Boolean IDs
    enum class BooleanID {NONE, MVAISOWP80, MVAISOWP90, MVAISOWPL, MVANOISOWP80, MVANOISOWP90, MVANOISOWPL, CUTBASEDHEEP};
    void SetBIDBit(BooleanID id, bool idbit);
    inline bool isMVAIsoWP80() const { return j_mvaIso_WP80; }
    inline bool isMVAIsoWP90() const { return j_mvaIso_WP90; }
    inline bool isMVAIsoWPLoose() const { return j_mvaIso_WPL; }
    inline bool isMVANoIsoWP80() const { return j_mvaNoIso_WP80; }
    inline bool isMVANoIsoWP90() const { return j_mvaNoIso_WP90; }
    inline bool isMVANoIsoWPLoose() const { return j_mvaNoIso_WPL; }
    inline bool isCutBasedHEEP() const { return j_cutBased_HEEP; }


    // cut-based ID
    enum class CutBasedID {NONE, CUTBASED};
    enum class WORKINGPOINT {NONE, VETO, LOOSE, MEDIUM, TIGHT};
    void SetCBIDBit(CutBasedID id, unsigned int value);
    inline WORKINGPOINT CutBased() const {return (WORKINGPOINT)j_cutBased; }

    // Private IDs
    bool Pass_CaloIdL_TrackIdL_IsoVL() const;
    bool Pass_HcToWABaseline() const;
    bool Pass_HcToWALooseRun2() const;
    bool Pass_HcToWALooseRun3() const;
    bool Pass_HcToWATight() const;
    // MVA scores
    enum class MVATYPE {NONE, MVAISO, MVANOISO, MVATTH};
    void SetMVA(MVATYPE type, float score);
    inline float MvaIso() const { return j_mvaIso; }
    inline float MvaNoIso() const { return j_mvaNoIso; }
    inline float MvaTTH() const { return j_mvaTTH; }

    // ID helper functions
    bool PassID(const TString ID) const;
    bool PassID(ElectronID ID) const;

private:
    // uncertainties
    float j_energyErr; 
    float j_dEsigmaUp, j_dEsigmaDown;
    // ID variables
    bool j_convVeto;
    unsigned char j_lostHits, j_seedGain, j_tightCharge;
    float j_sieie, j_hoe, j_eInvMinusPInv;
    float j_dr03EcalRecHitSumEt, j_dr03HcalDepth1TowerSumEt, j_dr03TkSumPt, j_dr03TkSumPtHEEP; // no j_e2x5OverE5x5, j_e1x5OverE5x5, j_calPFClusterIso, j_hcalPFClusterIso, j_dr03HcalTowerSumEt

    // SuperCluster
    float j_scEta, j_scPhi;
    float j_deltaEtaInSC, j_deltaEtaInSeed, j_deltaPhiInSC, j_deltaPhiInSeed; // missing phi and E?
    float j_ecalPFClusterIso, j_hcalPFClusterIso;

    // IDs
    bool j_mvaIso_WP80, j_mvaIso_WP90, j_mvaIso_WPL, j_mvaNoIso_WP80, j_mvaNoIso_WP90, j_mvaNoIso_WPL, j_cutBased_HEEP;
    unsigned char j_cutBased;
    
    float j_mvaIso, j_mvaNoIso, j_mvaTTH;

    // others
    float j_r9;
    float j_rho;
    short j_genPartIdx;
    unsigned char j_genPartFlav;
    short j_jetIdx;
    ClassDef(Electron, 1);
};

#endif

