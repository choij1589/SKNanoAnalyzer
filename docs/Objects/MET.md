# MET (Missing Transverse Energy)

The `METv` object represents Missing Transverse Energy with Type-I corrections applied using the latest JES/JER corrections.

## Contents

- [Overview](#overview)
- [Getting MET](#getting-met)
- [Type-I Correction](#type-i-correction)
- [Systematic Variations](#systematic-variations)
- [Usage Example](#usage-example)

---

## Overview

| Run   | Base MET Source | Algorithm                        |
| ----- | --------------- | -------------------------------- |
| Run 2 | CHS or PUPPI    | Selectable via `Event::MET_Type` |
| Run 3 | PUPPI           | Default                          |

---

## Getting MET

### From Event

```cpp
// Get raw MET vector
Particle MET_raw = ev.GetMETVector(Event::MET_Type::PUPPI, Event::MET_Syst::CENTRAL);

// Available MET types
enum class MET_Type { CHS, PUPPI };

// Available systematics
enum class MET_Syst { CENTRAL, UP, DOWN };
```

### With Type-I Correction

```cpp
// Apply Type-I correction (recommended)
Particle METv = ApplyTypeICorrection(MET_raw, allJets, allElectrons, allMuons,
                                      MyCorrection::variation::nom);
```

---

## Type-I Correction

The Type-I MET correction accounts for differences between raw and corrected object momenta.

### Method Signature

```cpp
Particle ApplyTypeICorrection(const Particle& MET,
                              const RVec<Jet>& jets,
                              const RVec<Electron>& electrons,
                              const RVec<Muon>& muons,
                              const MyCorrection::variation& unclustered_syst);
```

### Object Contributions

| Object    | Selection           | Correction Applied                    |
| --------- | ------------------- | ------------------------------------- |
| Jets      | pT > 15, \|η\| < 5  | Corrected − Raw momentum              |
| Muons     | pT > 5, \|η\| < 2.4 | Corrected − Raw (Rochester)           |
| Electrons | pT > 7, \|η\| < 2.5 | Corrected − Nominal (scale/smear)     |

### Unclustered Energy

- **PUPPI MET:** Uses NanoAOD stored `PuppiMET_pt/phiUnclusteredUp/Down`
- **CHS MET:** Calculated as MET + jets + leptons, with ±10% systematic

---

## Systematic Variations

### Correlated Variations

Pass scaled/smeared objects to propagate systematics:

```cpp
// JES/JER variation: pass varied jets
RVec<Jet> jets_jesUp = ...;  // Apply JES up variation
Particle METv_jesUp = ApplyTypeICorrection(MET_raw, jets_jesUp, electrons, muons,
                                            MyCorrection::variation::nom);

// Electron energy variation: pass varied electrons
RVec<Electron> electrons_up = ScaleElectrons(ev, electrons, "up");
Particle METv_eleUp = ApplyTypeICorrection(MET_raw, jets, electrons_up, muons,
                                            MyCorrection::variation::nom);

// Muon scale variation: pass varied muons
RVec<Muon> muons_up = ScaleMuons(muons, "up");
Particle METv_muUp = ApplyTypeICorrection(MET_raw, jets, electrons, muons_up,
                                           MyCorrection::variation::nom);
```

### Unclustered Energy Variation

```cpp
// Unclustered energy up/down
Particle METv_unclUp = ApplyTypeICorrection(MET_raw, jets, electrons, muons,
                                             MyCorrection::variation::up);
Particle METv_unclDn = ApplyTypeICorrection(MET_raw, jets, electrons, muons,
                                             MyCorrection::variation::down);
```

### XY Correction (Optional)

```cpp
// Apply XY correction on top of Type-I (not recommended for PUPPI)
myCorr->METXYCorrection(METv, ev.run(), ev.nPV(),
                        MyCorrection::XYCorrection_MetType::Type1PFMET);
```

---

## Usage Example

### DiLepton Analyzer Pattern

```cpp
void DiLepton::defineObjects() {
    // Get base MET
    Particle MET_default = ev.GetMETVector(Event::MET_Type::PUPPI, Event::MET_Syst::CENTRAL);

    // Get all objects (with corrections already applied)
    RVec<Jet> allJets = GetAllJets();
    RVec<Electron> allElectrons = GetAllElectrons();
    RVec<Muon> allMuons = GetAllMuons();

    // Apply Type-I correction
    Particle METv = ApplyTypeICorrection(MET_default, allJets, allElectrons, allMuons,
                                          MyCorrection::variation::nom);

    // Use METv for analysis
    float met_pt = METv.Pt();
    float met_phi = METv.Phi();
}
```

### Accessing Raw Object pT

The correction uses `OriginalPt()` method for raw momentum:

```cpp
// Jets
float jet_rawPt = jet.GetOriginalPt();

// Muons
float muon_rawPt = muon.OriginalPt();

// Electrons: Pt() is corrected in NanoAOD (raw = corrected for nominal)
```

---

## Notes

- Type-I correction removes dependency on NanoAOD's built-in correction
- Latest JES/JER corrections are always used, ensuring consistency with jet selections
- MET and object variations are 100% correlated when using scaled objects
- For systematic studies, see `DiLeptonSystematic.yaml` for proper variation names
