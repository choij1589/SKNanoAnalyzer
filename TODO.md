# TODO: PromptTreeProducer and MatrixTreeProducer Implementation

This document tracks the differences between the old Python-based SKFlatAnalyzer and the new C++ SKNanoAnalyzer implementations, highlighting missing features that need to be implemented.

## Completed Features ✅

### Core Framework
- [x] C++ implementation inheriting from TriLeptonBase
- [x] SystematicHelper integration for automatic systematic management
- [x] Separate TTree per systematic variation (matching old Python structure)
- [x] Support for both 1E2Mu and 3Mu channels
- [x] Event selection and kinematic cuts
- [x] MET Type-I correction application

### Systematic Variations
- [x] Weight-only systematics (using Central objects)
  - L1Prefire
  - PileupReweight
  - MuonIDSF
  - ElectronIDSF
  - EMuTrigSF (1E2Mu channel)
  - DblMuTrigSF (3Mu channel)
  - PileupJetIDSF (Run 2)
  - B-tagging systematics (HFcorr, HFuncorr, LFcorr, LFuncorr)
  - WZNjetsSF (Run 3 WZ/ZZ samples)

- [x] Object-variation systematics (requiring event re-processing)
  - JetEn (Up/Down)
  - JetRes (Up/Down)
  - ElectronEn (Up/Down)
  - ElectronRes (Up/Down)
  - MuonEn (Up/Down)
  - UnclusteredEn (Up/Down)

### Physics Features
- [x] Conversion sample filtering (DYJets, TTG, WWG)
  - Internal/external conversion lepton identification
  - Hadronic fake lepton rejection
- [x] Sample patching (DYJets vs ZGToLLG/DYGTo2LG)
  - Low pT region for DYJets
  - High pT region for ZGToLLG
- [x] Muon pair mass calculation (mass1, mass2)
- [x] Fake weight calculation (MatrixTreeProducer)
- [x] Fold calculation using Central METv (consistent across systematics)

### Scale Factors and Corrections
- [x] Muon reconstruction and ID scale factors
- [x] Electron reconstruction and ID scale factors
- [x] Trigger scale factors (EMu and DblMu)
- [x] B-tagging scale factors (Method 1a)
- [x] Pileup reweighting
- [x] L1 prefire weights
- [x] WZ Njets reweighting (Run 3)

### Tree Structure
- [x] Multiple TTrees: `Events_Central`, `Events_{Systematic}_Up`, `Events_{Systematic}_Down`
- [x] Branches: run, event, lumi, mass1, mass2, fold, weight
- [x] GraphNet score branches (12 branches: 3 signals × 4 classes)

## Theory Uncertainties ✅ COMPLETED

### Implementation Summary

Theory uncertainties are now implemented in `PromptTreeProducer` using pre-computed NanoAOD weights directly (no LHAPDF recalculation required).

**Enabled via:** `--userflags RunTheoryUnc`

#### PDF Uncertainties (100 variations) ✅
- [x] Create TTrees: `Events_PDF_0` through `Events_PDF_99`
- [x] Use NanoAOD `LHEPdfWeight[1-100]` directly

#### AlphaS Uncertainties (2 variations) ✅
- [x] Create TTrees: `Events_AlphaS_Up`, `Events_AlphaS_Down`
- [x] Use NanoAOD `LHEPdfWeight[101-102]` directly
- [ ] AlphaS factorization (requires LHAPDF) - placeholder for future

#### Scale Uncertainties (7 variations) ✅
- [x] Create TTrees: `Events_Scale_0,1,2,3,4,6,8` (skip 5,7)
- [x] Use NanoAOD `LHEScaleWeight[i]` directly

#### Parton Shower Uncertainties (4 variations) ✅
- [x] Create TTrees: `Events_PS_ISRUp`, `Events_PS_FSRUp`, `Events_PS_ISRDown`, `Events_PS_FSRDown`
- [x] Use NanoAOD `PSWeight[0-3]` directly

#### Implementation Details
- All theory uncertainties only run when `RunTheoryUnc` flag is set
- All use Central objects (weight-only variations)
- Skip theory trees for samples without weights (`HasTheoryWeights()` check)
- Total: 113 additional trees (100 PDF + 7 Scale + 4 PS + 2 AlphaS)

---

### 2. GraphNet Score Evaluation ✅ COMPLETED

**Implementation Summary:**
- [x] Load ParticleNet models in `initializeAnalyzer()` via `loadGraphNetModels()`
- [x] Implement `calculateFold()` function using Central METv for consistency
- [x] Implement `evalGraphNetScores()` function with full 4-class probability output
  - [x] Construct node features from muons, electrons, jets, bjets, METv
  - [x] Build k-NN graph (k=4) in C++ using brute force algorithm
  - [x] Run model inference for all 3 signal points
  - [x] Store all 4 class probabilities (signal, nonprompt, diboson, ttZ)
- [x] Integration with libtorch (PyTorch C++ API)
- [x] Model paths from SKNANO_DATA environment variable

**C++ Implementation:**
- TorchScript model loading via `torch::jit::load()`
- Graph input: 9 node features (E, Px, Py, Pz, Charge, IsMuon, IsElectron, IsJet, IsBjet)
- Edge construction: k-NN graph with k=4 computed in C++
- Era encoding: One-hot vector for Run2 eras (2016preVFP, 2016postVFP, 2017, 2018)
- Output: 12 branches per analyzer (3 signals × 4 classes)

**Key Fix Applied:**
- Models re-exported with CPPCompatibleParticleNet wrapper to eliminate torch_cluster dependency
- All DynamicEdgeConv layers now use pre-computed edge_index from C++
- No external Python dependencies required for inference

**Signal points:** MHc160_MA85, MHc130_MA90, MHc100_MA95
**Background classes:** nonprompt, diboson, ttZ

### 3. Code Organization Improvements

- [ ] Consider moving theory uncertainty handling to a separate helper class
- [ ] Add configuration file for model paths
- [ ] Add documentation for each systematic variation
- [ ] Add unit tests for systematic handling

### 4. Performance Optimization

- [ ] Profile systematic processing performance
- [ ] Consider caching Central objects more efficiently
- [ ] Optimize tree filling for many systematics

## Implementation Priority

1. **URGENT**: GraphNet score evaluation (required for analysis)
2. **HIGH**: Theory uncertainties (PDF, Scale, AlphaS, PS)
3. **MEDIUM**: Code organization improvements
4. **LOW**: Performance optimization

## Testing Checklist

When implementing missing features, verify:

- [ ] All TTrees are created correctly
- [ ] Fold is identical across all systematics
- [ ] GraphNet scores vary correctly with systematic variations
- [ ] Theory uncertainty weights are applied correctly
- [ ] Memory usage is acceptable with O(100) TTrees
- [ ] Output file structure matches old Python implementation

## Notes

- The old Python implementation is in:
  - `SKFlatAnalyzer/PyAnalyzers/PromptSkimmer.py`
  - `SKFlatAnalyzer/PyAnalyzers/MatrixSkimmer.py`

- The new C++ implementation is in:
  - `Analyzers/include/PromptTreeProducer.h`
  - `Analyzers/src/PromptTreeProducer.cc`
  - `Analyzers/include/MatrixTreeProducer.h`
  - `Analyzers/src/MatrixTreeProducer.cc`

- Both analyzers inherit from `TriLeptonBase`
- Systematics are managed by `SystematicHelper` reading from `TriLeptonSystematics.yaml`
