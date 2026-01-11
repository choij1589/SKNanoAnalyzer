# CLAUDE.md

## Project Overview

SKNanoAnalyzer is a C++ physics analysis framework for CMS experiment NanoAOD data. Supports Run 2 (2016-2018) and Run 3 (2022-2023) periods for dilepton studies, ttbar reconstruction, and systematic uncertainty evaluations.

## Quick Reference

```bash
source setup.sh              # Required once per session
./scripts/build.sh           # Full clean build
./scripts/rebuild.sh         # Incremental rebuild

# Job submission
SKNano.py -a AnalyzerName -i 'Sample*' -e 2022 -n 10
SKNano.py -a Skim_DiLepton -i DYJets -e 2022 --skimming_mode
```

## Directory Structure

```
SKNanoAnalyzer/
├── Analyzers/           # Analysis implementations (inherit AnalyzerCore)
├── AnalyzerTools/       # Utilities, corrections, systematic handling
├── DataFormats/         # Physics object classes
├── PyAnalyzers/         # Python postprocessing tools
├── python/              # Job submission (SKNano.py, sampleManager.py)
├── scripts/             # Build and utility scripts
├── external/            # Git submodules (RoccoR)
├── data/                # Era-specific corrections, efficiencies, samples
├── config/              # User config files (config.$USER)
├── SampleLists/         # Sample lists per era/analysis
├── templates/           # HTCondor job templates
└── CMakeLists.txt       # Root build configuration
```

## Class Hierarchy

**Physics Objects:**
```
Particle (base)
├── Lepton → Muon, Electron
├── Jet, FatJet
├── Tau, Photon
└── Gen*, LHE (generator-level)
```

**Analyzer Chain:**
```
SKNanoLoader → AnalyzerCore → [DiLeptonBase, VcbBase] → User Analyzers
```

## Creating Analyzers

1. Create `Analyzers/include/MyAnalyzer.h` inheriting from `AnalyzerCore`
2. Implement in `Analyzers/src/MyAnalyzer.cc`:
   ```cpp
   void MyAnalyzer::initializeAnalyzer() { /* setup */ }
   void MyAnalyzer::executeEvent() { /* per-event logic */ }
   ```
3. Add to `Analyzers/include/AnalyzersLinkDef.hpp`
4. Rebuild with `./scripts/rebuild.sh`

## Code Patterns

**Enum Classes (always prefer over strings/magic numbers):**
```cpp
enum class Muon::BooleanID { NONE, LOOSE, MEDIUM, TIGHT };
enum class Jet::JetID { TIGHT, LOOSE };
enum class MyCorrection::variation { nom, up, down };
```

**Object Selection:**
```cpp
RVec<Lepton*> leptons;
for (auto& mu : *muons) {
    if (mu.PassID(Muon::BooleanID::TIGHT) && mu.Pt() > 20.)
        leptons.push_back(&mu);
}
```

**Corrections:**
```cpp
MyCorrection::variation syst = variation::nom;  // or up/down
float sf = corrections.GetMuonIDSF(key, muon, syst);
```

## Systematic Uncertainties

Configured via YAML files in `AnalyzerTools/`:
- `DiLeptonSystematic.yaml`: Full systematic list
- `noSyst.yaml`: Nominal only (quick tests)

**Types:**
- **Event-loop**: Rerun event processing (JetRes, ElectronEn)
- **Weight-only**: Reweight without reprocessing (L1Prefire, TrigSF)

Use `SystematicHelper` class with `variation` enum (nom/up/down).

## Configuration Hierarchy

1. **User config**: `config/config.$USER` - package manager, paths, Telegram
2. **Sample config**: `data/{era}/Sample/CommonSampleInfo.json` - xsec, lumi
3. **Systematic config**: YAML files in `AnalyzerTools/`

## Build System

**CMake with optional Ninja:**
```bash
./scripts/build.sh           # Uses Make by default
./scripts/build.sh ninja     # Use Ninja (faster)
```

**Key Environment Variables** (set by setup.sh):
- `SKNANO_VERSION`: Data version (Run3_v13_Run2_v9)
- `SKNANO_BUILDDIR`, `SKNANO_INSTALLDIR`, `SKNANO_LIB`
- `LHAPDF_*`, `LIBTORCH_*`, `CORRECTION_*`, `ONNXRUNTIME_*`

**Dependencies** (conda/mamba):
- ROOT 6.32+, correctionlib, onnxruntime-cpp, PyTorch/libtorch
- LHAPDF (auto-installed via `scripts/install_lhapdf.sh`)

## Job Submission (HTCondor)

**DAG Workflow:** Analysis jobs → Hadd (merge) → Summary

```bash
# Basic submission
SKNano.py -a DiLepton -i 'DYJets*' -e 2022 -n 50

# With data reduction (process 1/N events)
SKNano.py -a DiLepton -i 'DYJets*' -e 2022 -n 50 --reduction 10

# Skimming mode (recommended before full analysis)
SKNano.py -a Skim_DiLepton -i DYJets -e 2022 --skimming_mode
```

**Options:**
- `-a`: Analyzer name
- `-i`: Sample pattern (regex supported)
- `-e`: Era (2022, 2022EE, 2023, 2023BPix, 2016preVFP, etc.)
- `-n`: Files per job
- `--skimming_mode`: Enable skim output
- `--reduction N`: Process 1/N events

## Multi-Era Support

| Run | Eras | NanoAOD |
|-----|------|---------|
| Run 2 | 2016preVFP, 2016postVFP, 2017, 2018 | v9 |
| Run 3 | 2022, 2022EE, 2023, 2023BPix        | v13 |

Era-specific data loaded automatically from `data/{SKNANO_VERSION}/{era}/`.

## CI/CD

GitHub Actions workflow (`.github/workflows/main.yml`):
- Builds on push/PR using Docker container
- Requires SSH key for RoccoR submodule (CERN GitLab)

## Common Tasks

**Add new physics object branch:**
1. Update `DataFormats/include/` class definition
2. Add branch in `AnalyzerTools/src/SKNanoLoader.cc`
3. Update `python/NanoAODv*.json` if needed

**Add new correction:**
1. Place correction file in `data/{era}/`
2. Implement getter in `AnalyzerTools/src/MyCorrection.cc`
3. Add systematic variation support if applicable

**Debug locally:**
```bash
root -l -b -q 'runAnalyzer.C("AnalyzerName", "input.root", "output.root")'
```
