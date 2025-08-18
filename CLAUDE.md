# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SKNanoAnalyzer is a C++ physics analysis framework for CMS (Compact Muon Solenoid) experiment data analysis. It processes NanoAOD files for Run 2 (2016-2018) and Run 3 (2022-2023) data periods, providing tools for dilepton studies, ttbar reconstruction, and systematic uncertainty evaluations.

## Build Commands

```bash
# Initial setup (required once per session)
source setup.sh

# Build entire project
./scripts/build.sh

# Incremental rebuild
./scripts/rebuild.sh
```

## Running Analysis

```bash
# Submit analysis jobs
SKNano.py -a AnalyzerName -i 'SamplePattern*' -e 2022 -n 10

# Example: Run DiLepton analyzer on DYJets samples
SKNano.py -a DiLepton -i 'DYJets*' -e 2022 -n 50 --reduction 10

# Skimming mode for data reduction
SKNano.py -a Skim_AnalyzerName -i SampleName -e 2022 --skimming_mode
```

## Architecture Overview

### Core Components
- **DataFormats/**: Physics object classes (Electron, Muon, Jet, etc.) with analysis methods
- **AnalyzerTools/**: Analysis utilities, corrections, systematic helpers, ML integration  
- **Analyzers/**: Specific analysis implementations inheriting from AnalyzerCore
- **External dependencies**: KinematicFitter, RoccoR, jsonpog-integration, LHAPDF, libtorch

### Key Classes
- **AnalyzerCore**: Base class providing common functionality for all analyzers
- **SKNanoLoader**: Handles NanoAOD file loading and branch management
- **SystematicHelper**: Manages systematic uncertainty variations
- **MyCorrection**: Centralized correction and scale factor application

### Data Organization
- Era-based structure: `data/{era}/` contains corrections, efficiencies, and sample metadata
- Sample information in `data/{era}/Sample/CommonSampleInfo.json`
- Systematic configurations in YAML files (DiLeptonSystematic.yaml, VcbSystematic.yaml)

## Development Patterns

### Creating New Analyzers
1. Create header in `Analyzers/include/` inheriting from `AnalyzerCore`
2. Implement in `Analyzers/src/` with required methods:
   - `initializeAnalyzer()`
   - `executeEvent()`
   - Virtual destructor
3. Add to `Analyzers/include/AnalyzersLinkDef.hpp`
4. Use enum classes for constants instead of magic numbers/strings

### Analyzer Implementation Structure
When implementing new analyzers, follow the established pattern from `DiLepton.cc`:

1. **Core Structure**:
   - `initializeAnalyzer()`: Set up channels, systematic helper, and analysis configuration
   - `executeEvent()`: Main event processing with sequential filters and systematic handling
   - Helper methods: `defineObjects()`, `selectEvent()`, `getWeights()`, `fillObjects()`

2. **Systematic Handling Pattern**:
   - Use `SystematicHelper` for managing systematic variations
   - Process Central objects first, then weight-only systematics
   - Handle object-variation systematics with event re-processing
   - Implement systematic-aware histogram organization: `channel/systematic/object/variable`

3. **Event Processing Flow**:
   - Apply sequential event filters (noise filter, veto maps, etc.)
   - Define physics objects with systematic variations
   - Apply channel-specific event selection
   - Calculate systematic weights and scale factors
   - Fill histograms with proper systematic categorization

4. **Object Definition Best Practices**:
   - Create object copies for systematic variations
   - Apply scale/smearing based on systematic names
   - Implement Type-I MET corrections
   - Sort objects by pT and apply selection criteria consistently

5. **Weight Calculation**:
   - Handle all systematic weight variations through `MyCorrection` class
   - Include trigger, lepton ID/reconstruction, b-tagging, and pileup scale factors
   - Apply generator weights, L1 prefire, and physics-specific reweighting

### Systematic Uncertainties
- Use `SystematicHelper` class for systematic variations
- Implement variations using `variation` enum (nom, up, down)
- Apply corrections through `MyCorrection` class methods

### Physics Object Usage
```cpp
// Example pattern for lepton selection
RVec<Lepton *> leptons;
TString MuonID = "Tight"; // Example ID
for (auto &mu : *muons) {
    if (mu.Pass(MuonID)) {
        leptons.push_back(&mu);
    }
}
```

## Dependencies and Environment

### Required Packages (managed via conda/mamba)
- ROOT 6.32+
- correctionlib (CMS corrections)
- onnxruntime-cpp (ML inference)
- PyTorch/libtorch (deep learning)
- LHAPDF (auto-installed via scripts/install\_lhapdf.sh)

### Configuration
- Personal config: `config/config.$USER` (system type, package manager, notifications)
- Sample lists: `SampleLists/{era}/` for different analysis categories

## Job Management

### HTCondor Integration
- Jobs submitted via HTCondor with DAG workflows
- Configurable memory (default 2GB) and CPU resources
- Automatic output file management and merging
- Telegram bot notifications (optional)

### Sample Selection
- Regex patterns for flexible sample matching
- Era-specific sample organization (2016preVFP, 2017, 2018, 2022, 2022EE, 2023, etc.)
- Cross-section and luminosity handled automatically

## Multi-Era Support

The framework supports unified analysis across data-taking periods:
- **Run 2**: 2016preVFP, 2016postVFP, 2017, 2018
- **Run 3**: 2022, 2022EE, 2023, 2023BPix

Era-specific corrections and configurations are automatically loaded based on the `-e` flag in job submission.
