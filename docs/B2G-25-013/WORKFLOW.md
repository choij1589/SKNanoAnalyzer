# Workflow for B2G-25-013

## 1. Object Reviews (DiLepton)
Run `DiLepton` analyzer to check object kinematics and basic distributions.
- **Script**: `scripts/submission/DiLepton.sh`
- **Channels**: `RunDiMu`, `RunEMu`
- **Flags**: `RunSyst`
- **Inputs**:
  - Data Streams: `MuonEG` (for EMu), `DoubleMuon,Muon,Muon0,Muon1` (for DiMu)
  - MC: `DYJets`, `WJets`, `TTLL_powheg`, `TTLJ_powheg`
  - Sample Lists: `DiLepton.txt` (Run2NanoV9 or Run3NanoV13)

## 2. Background Estimation
### Fake Rate Measurement
- `MeasFakeRateV3`
- Closure test: `ClosFakeRate`

## 3. Control Regions (Conversion Rate and WZNj Reweighting)
Run `PromptAnalyzer` and `MatrixAnalyzer` for background estimation and reweighting derivations.
- **Script**: `scripts/submission/TriLepton.sh`
- **Mode**: `CR`
- **Channels**: `Run1E2Mu`, `Run2E1Mu`, `Run3Mu`
- **Flags**: `RunCR`, `NoTreeMode`, `RunSyst`
- **Inputs**: Data, `TriLepton.txt`, `SignalSamples.txt`

### WZ Normalization
Check WZ scaling factors.
- **Script**: `scripts/submission/TriLepton.sh`
- **Mode**: `WZ`
- **Target**: `WZTo3LNu` sample
- **Comparison**: `RunSyst` vs `RunNoWZSF`

## 4. Signal Estimation (Signal Region)
Run final analysis on Signal Region.
- **Script**: `scripts/submission/TriLepton.sh`
- **Mode**: `SR`
- **Channels**: `Run1E2Mu`, `Run2E1Mu`, `Run3Mu`
- **Flags**: `RunSyst`, `RunTheoryUnc` (for signals)
- **Inputs**:
  - Data: `Skim_TriLep` datasets
  - MC: `TriLepton.txt`, `SignalSamples.txt`
