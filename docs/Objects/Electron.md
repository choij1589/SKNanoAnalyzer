# Electron

The `Electron` class represents reconstructed electrons. Inherits from `Lepton` → `Particle` → `TLorentzVector`.

## Contents

- [Class Hierarchy](#class-hierarchy)
- [Getting Electrons](#getting-electrons)
- [Identification](#identification)
- [Isolation](#isolation)
- [Energy Corrections](#energy-corrections)
- [Scale Factors](#scale-factors)
- [Properties Reference](#properties-reference)
- [Run 2 vs Run 3](#run-2-vs-run-3)

---

## Class Hierarchy

```text
TLorentzVector
└── Particle        (charge, 4-momentum)
    └── Lepton      (dXY, dZ, IP3D, isolation)
        └── Electron
```

**Files:** `DataFormats/include/Electron.h`, `DataFormats/src/Electron.cc`

---

## Getting Electrons

### From AnalyzerCore

```cpp
// All electrons (GAP region excluded automatically)
RVec<Electron> allElectrons = GetAllElectrons();

// With ID, pT, eta cuts
RVec<Electron> electrons = GetElectrons("POGTight", 20., 2.5);

// With HEM veto (2018 only)
RVec<Electron> electrons = GetElectrons("POGTight", 20., 2.5, true);

// Select from existing collection (enum-based)
RVec<Electron> selected = SelectElectrons(electrons, Electron::ElectronID::POG_TIGHT, 25., 2.4);
```

---

## Identification

### Enum-Based IDs (Recommended)

```cpp
enum class ElectronID {
    NOCUT,
    POG_VETO, POG_LOOSE, POG_MEDIUM, POG_TIGHT,
    POG_HEEP,
    POG_MVAISO_WP80, POG_MVAISO_WP90, POG_MVAISO_WPL,
    POG_MVANOISO_WP80, POG_MVANOISO_WP90, POG_MVANOISO_WPL,
};

if (electron.PassID(Electron::ElectronID::POG_TIGHT)) { ... }
```

### String-Based IDs

| String              | Description                       |
| ------------------- | --------------------------------- |
| `"POGVeto"`         | Cut-based Veto                    |
| `"POGLoose"`        | Cut-based Loose                   |
| `"POGMedium"`       | Cut-based Medium                  |
| `"POGTight"`        | Cut-based Tight                   |
| `"POGHEEP"`         | High-Energy Electron Pairs        |
| `"POGMVAIsoWP80"`   | MVA with Iso, 80% eff WP          |
| `"POGMVAIsoWP90"`   | MVA with Iso, 90% eff WP          |
| `"POGMVANoIsoWP80"` | MVA without Iso, 80% eff WP       |
| `"POGMVANoIsoWP90"` | MVA without Iso, 90% eff WP       |
| `"HcToWATight"`     | Analysis-specific tight           |
| `"HcToWALooseRun2"` | Analysis-specific loose (Run 2)   |
| `"HcToWALooseRun3"` | Analysis-specific loose (Run 3)   |

### MVA Scores

```cpp
float mvaIso   = electron.MvaIso();       // MVA with isolation
float mvaNoIso = electron.MvaNoIso();     // MVA without isolation
float mvaTTH   = electron.MvaTTH();       // ttH-optimized

// Boolean WP checks
bool pass = electron.isMVAIsoWP80();
bool pass = electron.isMVANoIsoWP90();
```

### Cut-Based Working Points

```cpp
enum class WORKINGPOINT { NONE, VETO, LOOSE, MEDIUM, TIGHT };

// Check if passes at least Medium
if ((int)electron.CutBased() >= (int)Electron::WORKINGPOINT::MEDIUM) { ... }
```

---

## Isolation

Inherited from `Lepton`:

```cpp
float tkRelIso = electron.TkRelIso();        // Track-based
float pfIso03  = electron.PfRelIso03();      // PF ΔR < 0.3
float pfIso04  = electron.PfRelIso04();      // PF ΔR < 0.4
float miniIso  = electron.MiniPFRelIso();    // Mini-isolation (pT-dependent cone)
```

---

## Energy Corrections

### Scale and Smearing

Energy corrections are pre-applied in NanoAOD. Use these for systematic variations:

```cpp
// Scale variation (data-like uncertainty)
RVec<Electron> scaled = ScaleElectrons(event, electrons, "up");    // or "down"

// Smearing variation (MC resolution uncertainty)
RVec<Electron> smeared = SmearElectrons(electrons, "up");          // or "down", "nom" (Run 3)
```

### Energy Uncertainty

```cpp
float energyErr = electron.energyErr();

// Run 2: resolution variations stored in NanoAOD
float dEup   = electron.dEsigmaUp();
float dEdown = electron.dEsigmaDown();
```

### Direct Correction Access

```cpp
// Scale uncertainty
float scale = myCorr->GetElectronScaleUnc(
    electron.scEta(), electron.SeedGain(), run, electron.r9(), electron.Pt(),
    MyCorrection::variation::up);

// Smearing (Run 3)
float smear = myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::nom, seed);
```

---

## Scale Factors

```cpp
// Reconstruction SF
float recoSF = myCorr->GetElectronRECOSF(eta, pt, phi, MyCorrection::variation::nom);

// Identification SF
float idSF = myCorr->GetElectronIDSF("Tight", eta, pt, phi, MyCorrection::variation::nom);

// Trigger SF
float trigSF = myCorr->GetElectronTriggerSF("Ele32", eta, pt, phi, MyCorrection::variation::nom);
```

---

## Properties Reference

### Supercluster

```cpp
float scEta = electron.scEta();

// Eta regions
enum class ETAREGION { IB, OB, GAP, EC };  // |η|: <0.8, <1.444, <1.566, else
ETAREGION region = electron.etaRegion();

// Track-SC matching
float dEtaSC   = electron.deltaEtaInSC();
float dPhiSC   = electron.deltaPhiInSC();
float dEtaSeed = electron.deltaEtaInSeed();
float dPhiSeed = electron.deltaPhiInSeed();
```

### Shower Shape

```cpp
float sieie = electron.sieie();               // σ_iηiη
float hoe   = electron.hoe();                 // H/E
float r9    = electron.r9();                  // R9
float eInvMinusPInv = electron.eInvMinusPInv();
```

### Track Quality

```cpp
bool  convVeto    = electron.ConvVeto();      // Conversion veto
uchar lostHits    = electron.LostHits();      // Missing inner hits
uchar tightCharge = electron.TightCharge();   // Charge consistency (0/1/2)
uchar seedGain    = electron.SeedGain();
```

### Isolation Sums (ΔR < 0.3)

```cpp
float ecalSum   = electron.dr03EcalRecHitSumEt();
float hcalSum   = electron.dr03HcalDepth1TowerSumEt();
float tkSum     = electron.dr03TkSumPt();
float tkSumHEEP = electron.dr03TkSumPtHEEP();

// PF cluster isolation
float ecalPF = electron.ecalPFClusterIso();
float hcalPF = electron.hcalPFClusterIso();
```

### Impact Parameters (from Lepton)

```cpp
float dxy   = electron.dXY();
float dz    = electron.dZ();
float ip3d  = electron.IP3D();
float sip3d = electron.SIP3D();
```

### Generator Matching (MC)

```cpp
short genIdx  = electron.GenPartIdx();    // Index in GenPart
uchar genFlav = electron.GenPartFlav();   // Generator flavor
short jetIdx  = electron.JetIdx();        // Matched jet index
```

---

## Run 2 vs Run 3

| Aspect             | Run 2 (NanoAODv9)        | Run 3 (NanoAODv13)             |
| ------------------ | ------------------------ | ------------------------------ |
| ID Algorithm       | Fall17v2                 | Winter22v1                     |
| MVA branches       | `mvaFall17V2*`           | `mvaIso*`, `mvaNoIso*`         |
| Energy smearing    | Stored in NanoAOD        | Applied via correction         |
| φ in scale factors | Not used                 | Used in some SFs               |
| Eras               | 2016pre/post, 2017, 2018 | 2022, 2022EE, 2023, 2023BPix   |

> **Note:** NanoAODv12 has a bug in MVA NoIso ID. Use v13 for Run 3. See [MVAElectronIDRun3](https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun3).

---

## HEM Veto (2018)

The 2018 data has a HEM detector issue affecting electrons in η < -1.25 and -1.62 < φ < -0.82:

```cpp
// Check manually
bool isHEM = IsHEMElectron(electron);

// Or apply during selection
RVec<Electron> electrons = GetElectrons("POGTight", 20., 2.4, true);  // vetoHEM=true
```

---

## References

- [CMS EGamma POG Run2](https://cms-egamma.github.io/Run2/)
- [EGM Scale Factors JSON](https://cms-egamma.github.io/EgammaSFJSON/)
- [Electron Performance 2022-2023](https://cds.cern.ch/record/2904365)
- [correctionlib Framework](https://github.com/cms-nanoAOD/correctionlib)
