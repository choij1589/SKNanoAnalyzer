# METv (Missing Transverse Energy Vector)

## Overview

The `METv` object in SKNanoAnalyzer represents the Missing Transverse Energy corrected with Type-I corrections using the latest Jet Energy Scale (JES) and Jet Energy Resolution (JER) corrections. This ensures consistency between the MET calculation and the jet corrections applied in the analysis.

## Implementation

### Base MET Source
- **Run 3**: PUPPI MET (`Event::MET_Type::PUPPI`)
- **Run 2**: CHS MET / PUPPI MET (Can be selected by `Event::MET_type::CHS/PUPPI`)

### Type-I Correction Method

The Type-I MET correction is implemented in `AnalyzerCore::ApplyTypeICorrection()` and application example can be found in the `DiLepton::defineObjects()` method. The correction accounts for:

#### Object Corrections
1. **Jets** (pT > 15 GeV, |η| < 5.0 ):
   - Subtracts difference between corrected and raw jet momentum
   - Uses latest JES/JER corrections from framework

2. **Muons** (pT > 5 GeV, |η| < 2.4):
   - Subtracts difference between corrected and raw muon momentum
   - Uses RoccoR method

3. **Electrons** (pT > 7 GeV, |η| < 2.5):
   - Raw = corrected in NanoAOD for nominal case
   - Only systematic variations are applied on the fly

4. **Unclustered Energy**:
   - For PUPPI MET, use NanoAOD stored unclustered energy
   - For CHS MET, calculated as: MET + jets + leptons, ±10% systematic variation applied

### Systematic Variations

The MET correction uses scaled objects from `GetAll*` method. The scaled objects are propagated to METv when calling ApplyTypeICorrection,
so METv and object scale variations are 100% correlated.
For unclustered energy up / down, it should be done by varing the name of systematic, Check `DiLeptonSystematic.yaml`.


### Technical Implementation

#### Method Signature
```cpp
Particle ApplyTypeICorrection(const Particle& MET,
                              const RVec<Jet>& jets,
                              const RVec<Electron>& electrons,
                              const RVec<Muon>& muons,
                              const MyCorrection::variation& var);
```

#### Systematic Handling
- Correlated variations (JES/JER, MuonEn, ElectronEn/Res) are applied by passing the already scaled/smeared `jets`, `muons`, and `electrons`.
- For unclustered energy, the systematic variations are handled by using `PuppiMET_pt/phiUnclusterdUp/Down` in NanoAOD branches for PUPPI MET, or
  re-calculated from other objects for CHS MET.

#### Object Access
- **Jets**: Uses `OriginalPt()` method for raw momentum
- **Muons**: Uses `OriginalPt()` method (renamed from `MiniAODPt()`)
- **Electrons**: `Pt()`: corrected in NanoAOD

### Usage in Analysis

#### DiLepton Analyzer
```cpp
// In defineObjects method
Particle MET_default = ev.GetMETVector(Event::MET_Type::PUPPI, Event::MET_Syst::CENTRAL);
Particle METv = ApplyTypeICorrection(MET, allJets, allElectrons, allMuons);

// XY correction is applied on top of Type-I correction
// Not recommended for PUPPI MET
myCorr->METXYCorrection(METv, ev.run(), ev.nPV(), MyCorrection::XYCorrection_MetType::Type1PFMET);
```

#### Systematic Processing
The MET correction is automatically applied for each systematic variation:
- Central: Uses nominal object corrections
- Systematic variations: Uses appropriate object variations and sources

### Notes

- The correction removes the dependency on NanoAOD's built-in Type-I correction
- Latest corrections are always used, ensuring consistency with jet selections
- Separate handling for Run 2 vs Run 3 data (future enhancement)
- Unclustered energy systematic properly accounts for object contributions
