# Muon

The `Muon` class represents reconstructed muons. Inherits from `Lepton` → `Particle` → `TLorentzVector`.

## Contents

- [Class Hierarchy](#class-hierarchy)
- [Getting Muons](#getting-muons)
- [Identification](#identification)
- [Isolation](#isolation)
- [Momentum Corrections](#momentum-corrections)
- [Scale Factors](#scale-factors)
- [Properties Reference](#properties-reference)
- [Run 2 vs Run 3](#run-2-vs-run-3)

---

## Class Hierarchy

```text
TLorentzVector
└── Particle        (charge, 4-momentum)
    └── Lepton      (dXY, dZ, IP3D, isolation)
        └── Muon
```

**Files:** `DataFormats/include/Muon.h`, `DataFormats/src/Muon.cc`

---

## Getting Muons

### From AnalyzerCore

```cpp
// All muons with Rochester corrections applied
RVec<Muon> allMuons = GetAllMuons();

// With ID, pT, eta cuts
RVec<Muon> muons = GetMuons("POGTight", 20., 2.4);

// Select from existing collection (enum-based)
RVec<Muon> selected = SelectMuons(muons, Muon::MuonID::POG_TIGHT, 25., 2.4);
```

---

## Identification

### Enum-Based IDs (Recommended)

```cpp
enum class MuonID {
    NOCUT,
    // Boolean IDs
    POG_TIGHT, POG_MEDIUM, POG_MEDIUM_PROMPT, POG_LOOSE,
    POG_SOFT, POG_SOFT_MVA, POG_TRIGGER_LOOSE,
    // High-pT IDs
    POG_TRACKER_HIGH_PT, POG_GLOBAL_HIGH_PT,
    // Mini Isolation
    POG_MINISO_LOOSE, POG_MINISO_MEDIUM, POG_MINISO_TIGHT, POG_MINISO_VTIGHT,
    // Multi Isolation
    POG_MULTISO_LOOSE, POG_MULTISO_MEDIUM,
    // MVA Muon
    POG_MVA_MU_MEDIUM, POG_MVA_MU_TIGHT,
    // PF Isolation
    POG_PFISO_VLOOSE, POG_PFISO_LOOSE, POG_PFISO_MEDIUM,
    POG_PFISO_TIGHT, POG_PFISO_VTIGHT, POG_PFISO_VVTIGHT,
    // PUPPI Isolation
    POG_PUPPIISO_LOOSE, POG_PUPPIISO_MEDIUM, POG_PUPPIISO_TIGHT,
    // Tracker Isolation
    POG_TKISO_LOOSE, POG_TKISO_TIGHT,
};

if (muon.PassID(Muon::MuonID::POG_TIGHT)) { ... }
```

### String-Based IDs

| String               | Description                            |
| -------------------- | -------------------------------------- |
| `"POGTight"`         | Cut-based Tight                        |
| `"POGMedium"`        | Cut-based Medium                       |
| `"POGMediumPrompt"`  | Medium with prompt cuts                |
| `"POGLoose"`         | Cut-based Loose                        |
| `"POGSoft"`          | Soft muon ID                           |
| `"POGSoftMVA"`       | MVA-based soft ID                      |
| `"POGTriggerLoose"`  | Loose trigger ID                       |
| `"POGTrackerHighPt"` | High-pT tracker muon                   |
| `"POGGlobalHighPt"`  | High-pT global muon                    |
| `"POGMiniIso*"`      | Mini isolation (Loose/Medium/Tight/VTight) |
| `"POGMultiIso*"`     | Multi isolation (Loose/Medium)         |
| `"POGPfIso*"`        | PF isolation (VLoose to VVTight)       |
| `"POGPuppiIso*"`     | PUPPI isolation (Loose/Medium/Tight)   |
| `"POGTkIso*"`        | Tracker isolation (Loose/Tight)        |
| `"HcToWATight"`      | Analysis-specific tight                |
| `"HcToWALoose"`      | Analysis-specific loose                |

### MVA Scores

```cpp
float softMva  = muon.SoftMva();      // Soft MVA score
float mvaLowPt = muon.MvaLowPt();     // Low-pT MVA
float mvaTTH   = muon.MvaTTH();       // ttH-optimized
```

### Muon Type

```cpp
bool isTracker    = muon.isTracker();
bool isStandalone = muon.isStandalone();
bool isGlobal     = muon.isGlobal();
```

---

## Isolation

### Working Point IDs

```cpp
enum class WorkingPoint { NONE, VLOOSE, LOOSE, MEDIUM, TIGHT, VTIGHT, VVTIGHT };

// Check working point level
if ((int)muon.PfIsoId() >= (int)Muon::WorkingPoint::TIGHT) { ... }
if ((int)muon.MiniIsoId() >= (int)Muon::WorkingPoint::MEDIUM) { ... }
```

### Available Isolation Types

```cpp
WorkingPoint pfIso    = muon.PfIsoId();       // PF-based (6 WPs)
WorkingPoint miniIso  = muon.MiniIsoId();     // Mini isolation (4 WPs)
WorkingPoint multiIso = muon.MultiIsoId();    // Combined (2 WPs)
WorkingPoint puppiIso = muon.PuppiIsoId();    // PUPPI-based (3 WPs)
WorkingPoint tkIso    = muon.TkIsoId();       // Tracker-only (2 WPs)
WorkingPoint mvaMu    = muon.MvaMuId();       // MVA muon (2 WPs)
WorkingPoint highPt   = muon.HighPtId();      // High-pT (1=tracker, 2=global)
```

### Raw Isolation Values (from Lepton)

```cpp
float tkRelIso = muon.TkRelIso();        // Track-based
float pfIso03  = muon.PfRelIso03();      // PF ΔR < 0.3
float pfIso04  = muon.PfRelIso04();      // PF ΔR < 0.4
float miniIso  = muon.MiniPFRelIso();    // Mini-isolation
```

---

## Momentum Corrections

### Rochester Corrections

Rochester corrections are applied automatically in `GetAllMuons()`. For systematic variations:

```cpp
// Scale variation
RVec<Muon> scaled = ScaleMuons(muons, "up");    // or "down"
```

### Direct Correction Access

```cpp
// Data
float roccor = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::nom);

// MC (with gen-matching)
float matched_pt = matched_gen.Pt();
float roccor = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::nom, matched_pt);

// Systematic uncertainty
float roccor_err = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::up) - roccor;
```

### Momentum Properties

```cpp
float originalPt = muon.OriginalPt();           // Pre-Rochester pT
float scaleUp    = muon.MomentumScaleUp();      // Scale up variation
float scaleDown  = muon.MomentumScaleDown();    // Scale down variation
int   nLayers    = muon.nTrackerLayers();       // Tracker layers
```

---

## Scale Factors

```cpp
// Reconstruction SF (Run 2 only)
float recoSF = myCorr->GetMuonRECOSF(muon, MyCorrection::variation::nom);

// Identification SF
float idSF = myCorr->GetMuonIDSF("NUM_TightID_DEN_TrackerMuons", muon,
                                  MyCorrection::variation::nom);

// Isolation SF
float isoSF = myCorr->GetMuonIDSF("NUM_TightRelIso_DEN_TightIDandIPCut", muon,
                                   MyCorrection::variation::nom);

// Trigger SF
float trigSF = myCorr->GetMuonTriggerSF("IsoMu24", muons, MyCorrection::variation::nom);
```

---

## Properties Reference

### Impact Parameters (from Lepton)

```cpp
float dxy   = muon.dXY();
float dz    = muon.dZ();
float ip3d  = muon.IP3D();
float sip3d = muon.SIP3D();
```

### Generator Matching (MC)

```cpp
short genIdx  = muon.GenPartIdx();    // Index in GenPart
uchar genFlav = muon.GenPartFlav();   // Generator flavor
short jetIdx  = muon.JetIdx();        // Matched jet index
```

---

## Run 2 vs Run 3

| Aspect               | Run 2 (NanoAODv9)        | Run 3 (NanoAODv13)             |
| -------------------- | ------------------------ | ------------------------------ |
| Rochester Corr.      | RoccoR package           | RoccoR (preliminary)           |
| RECO SF              | Applied via POG SFs      | Not needed                     |
| Isolation            | PF/Mini/PUPPI            | Improved PUPPI-based           |
| Systematic Framework | Manual                   | correctionlib                  |
| Eras                 | 2016pre/post, 2017, 2018 | 2022, 2022EE, 2023, 2023BPix   |

> **Note:** Run 3 Rochester corrections are preliminary.

---

## References

- [CMS Muon POG](https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsMUO)
- [Muon POG Recommendations](https://muon-wiki.docs.cern.ch/guidelines/corrections/)
- [Rochester Corrections](https://gitlab.cern.ch/hyseo/RoccoR)
- [Run 3 Muon Performance](https://www.researchgate.net/publication/377547277)
