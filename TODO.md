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

## Missing Features ❌

### 1. Theory Uncertainties (HIGH PRIORITY)

Theory uncertainties are completely missing from the C++ implementation. The old Python version had full support.

#### PDF Uncertainties
- [ ] PDF reweighting with 100 member variations
- [ ] Create TTrees: `Events_PDFReweight_0` through `Events_PDFReweight_99`
- [ ] Integration with LHAPDF
- [ ] Use `pdfReweight` object from AnalyzerCore

**Old Python implementation:**
```python
self.NPDF = 100
for pdf_idx in range(self.NPDF):
    pdf_syst = f"PDFReweight_{pdf_idx}"
    self.weight[pdf_syst][0] = w_norm * ... * self.weight_PDF.at(pdf_idx)
    self.tree[pdf_syst].Fill()
```

**Required for C++:**
- Access PDF weights via `pdfReweight->GetPDFWeight(pdf_idx)`
- Create 100 additional TTrees
- Apply PDF weights to Central event (no object variation)

#### AlphaS Uncertainties
- [ ] AlphaS variations (2 members: down, up)
- [ ] AlphaS factorization variations (2 members: down, up)
- [ ] Create TTrees: `Events_AlpS_down`, `Events_AlpS_up`, `Events_AlpSfact_down`, `Events_AlpSfact_up`

**Old Python implementation:**
```python
self.alphas_variations = ["AlpS_down", "AlpS_up", "AlpSfact_down", "AlpSfact_up"]
self.weight["AlpS_down"][0] = w_norm * ... * self.weight_AlphaS.at(0)
```

**Required for C++:**
- Access AlphaS weights from LHE information
- Create 4 additional TTrees
- Apply weights to Central event

#### Scale Uncertainties
- [ ] Renormalization and factorization scale variations (9 variations, skip indices 5 and 7)
- [ ] Create TTrees: `Events_ScaleVar_0` through `Events_ScaleVar_8` (excluding 5 and 7)

**Old Python implementation:**
```python
self.NSCALE = 9
for scale_idx in range(self.NSCALE):
    if scale_idx == 5 or scale_idx == 7: continue
    scale_syst = f"ScaleVar_{scale_idx}"
    self.weight[scale_syst][0] = w_norm * ... * self.weight_Scale.at(scale_idx)
```

**Required for C++:**
- Access scale weights via `GetScaleVariation(muF_syst, muR_syst)`
- Create 7 additional TTrees (9 - 2 skipped)
- Apply weights to Central event

#### Parton Shower Uncertainties
- [ ] PS ISR/FSR variations (4 variations)
- [ ] Create TTrees: `Events_PSVar_0` through `Events_PSVar_3`

**Old Python implementation:**
```python
self.NPSSYST = 4
for ps_idx in range(self.NPSSYST):
    ps_syst = f"PSVar_{ps_idx}"
    self.weight[ps_syst][0] = w_norm * ... * self.weight_PSSyst.at(ps_idx)
```

**Required for C++:**
- Access PS weights via `GetPSWeight(ISR_syst, FSR_syst)`
- Create 4 additional TTrees
- Apply weights to Central event

#### Implementation Notes
- All theory uncertainties should only run when `RunTheoryUnc` flag is set
- All use Central objects (no object variations)
- All modify weights only
- Need to check if theory weights are available in the sample

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

## Missing filter efficiencies for signal samples
/TTToHcToWAToMuMu_MHc-70_MA-65_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM
/TTToHcToWAToMuMu_MHc-160_MA-50_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM
/TTToHcToWAToMuMu_MHc-115_MA-27_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v3/MINIAODSIM
/TTToHcToWAToMuMu_MHc-100_MA-95_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
/TTToHcToWAToMuMu_MHc-115_MA-87_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM
/TTToHcToWAToMuMu_MHc-85_MA-15_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM
