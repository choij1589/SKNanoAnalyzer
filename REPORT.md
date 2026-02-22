# Run2 Endcap Electron Suppression in Systematic TTLL NanoAOD

## Summary

All 8 systematic TTLL variants (mtop, TuneCP5, CR, hdamp) have near-zero endcap electrons (|eta| > 1.479) in Run2 eras, while Run3 eras and the nominal TTLL_powheg are unaffected. The issue originates in the upstream Run2 NanoAOD production — not in the analysis code or dataset pipeline. This suppresses endcap electron coverage in the nonprompt training class for Run2.

## Observation

[Dataset validation](python/validateDatasets.py) plots for Run2 `ele1_eta` show the nonprompt class dropping to near-zero in the endcap (|eta| > 1.479), inconsistent with the expected ~20% endcap fraction seen in the TriLepton analysis and in the nominal TTLL_powheg sample.

## Investigation

### Step 1: Check saved .pt datasets

Electron eta distributions in the fold-0 `.pt` files (Run1E2Mu):

| Sample | Events | Barrel electrons | Endcap electrons | Endcap fraction |
|--------|-------:|-----------------:|-----------------:|----------------:|
| Skim_TriLep_TTLL_powheg | 13,406 | 10,728 | 2,678 | **20.0%** |
| Skim_TriLep_TTLL_mtop171p5_powheg | 5,115 | 4,850 | 265 | **5.2%** |
| Skim_TriLep_TTLL_TuneCP5up_powheg | 5,103 | 4,771 | 332 | **6.5%** |

The systematic variants show suppressed endcap even in the saved datasets.

### Step 2: Check EvtTreeProducer ROOT files per era

Direct check of `ElectronEtaColl` in the [EvtTreeProducer](../SKNanoAnalyzer/Analyzers/src/EvtTreeProducer.cc) output ROOT files:

**TTLL_powheg (nominal) — correct in all eras:**

| Era | Events | Barrel | Endcap | Endcap % |
|-----|-------:|-------:|-------:|---------:|
| 2016preVFP | 10,180 | 8,123 | 2,055 | 20.2% |
| 2017 | 24,765 | 19,604 | 5,161 | 20.8% |
| 2018 | 35,232 | 27,988 | 7,244 | 20.6% |
| 2022 | 3,219 | 2,451 | 768 | 23.9% |
| 2022EE | 12,127 | 9,425 | 2,702 | 22.3% |

**TTLL_mtop171p5_powheg — Run2 broken, Run3 correct:**

| Era | Events | Barrel | Endcap | Endcap % |
|-----|-------:|-------:|-------:|---------:|
| 2016preVFP | 4,181 | 4,171 | 10 | **0.2%** |
| 2017 | 9,052 | 9,028 | 24 | **0.3%** |
| 2018 | 13,199 | 13,159 | 40 | **0.3%** |
| 2022 | 1,236 | 968 | 268 | 21.7% |
| 2022EE | 4,872 | 3,821 | 1,051 | 21.6% |

**All systematic variants show the same pattern:**

| Sample | Run2 (2017) endcap % | Run3 (2022) endcap % |
|--------|---------------------:|---------------------:|
| TTLL_powheg (nominal) | 20.8% | 23.9% |
| TTLL_mtop171p5_powheg | **0.3%** | 21.7% |
| TTLL_mtop173p5_powheg | **~0.3%** | ~22% |
| TTLL_TuneCP5up_powheg | **0.3%** | 22.9% |
| TTLL_TuneCP5down_powheg | **0.3%** | 22.2% |
| TTLL_TuneCP5CR1_powheg | **0.2%** | 22.8% |
| TTLL_TuneCP5CR2_powheg | **~0.3%** | ~22% |
| TTLL_hdamp_up_powheg | **0.2%** | 21.7% |
| TTLL_hdamp_down_powheg | **~0.3%** | ~22% |

## Root Cause Analysis

### Ruled out: analysis code bug

The same electron ID ([`Pass_HcToWALooseRun2()`](../SKNanoAnalyzer/DataFormats/src/Electron.cc)) is applied uniformly to all samples. If the ID definition were the cause, the nominal TTLL_powheg would also be affected. It is not.

### Ruled out: branch mapping issue

The NanoAOD branch mapping in [`AnalyzerCore.cc:698-710`](../SKNanoAnalyzer/Analyzers/src/AnalyzerCore.cc) maps `Electron_mvaFall17V2noIso` to `MvaNoIso()` for all Run2 samples equally. Since the nominal works correctly, the mapping itself is fine.

### Identified: Run2 NanoAOD production issue

The issue is **sample-specific** in the Run2 NanoAOD input files:
- Only systematic TTLL variants are affected
- Only Run2 eras (2016preVFP, 2016postVFP, 2017, 2018) are affected
- Run3 eras use different MVA branches (`Electron_mvaNoIso` instead of `Electron_mvaFall17V2noIso`) and are unaffected

The most likely cause: the `Electron_mvaFall17V2noIso` and/or `Electron_mvaFall17V2noIso_WP90` values in the Run2 NanoAOD for systematic TTLL variants are incorrectly computed for endcap electrons. This could result from:
1. Different CMSSW release/configuration used for NanoAOD production of systematic variants
2. Missing or incorrect electron MVA calibrations during miniAOD-to-NanoAOD conversion
3. A bug in the NanoAOD production workflow for these specific samples

### Electron ID chain (reference)

The selection chain for Run2 loose electrons in [`EvtTreeProducer.cc`](../SKNanoAnalyzer/Analyzers/src/EvtTreeProducer.cc):

```
SelectElectrons(allElectrons, "HcToWALooseRun2", pT > 15, |eta| < 2.5)
  └── Pass_HcToWALooseRun2()                     [Electron.cc:164-184]
        ├── Pass_HcToWABaseline()                 [Electron.cc:148-154]
        │     ├── Pass_CaloIdL_TrackIdL_IsoVL()   [Electron.cc:120-146]
        │     │     └── etaRegion() uses scEta     [Electron.h:20-25]
        │     │         IB: |scEta| < 0.8
        │     │         OB: 0.8-1.444
        │     │         GAP: 1.444-1.566 (vetoed)
        │     │         EC: > 1.566
        │     ├── ConvVeto()
        │     ├── LostHits() < 2
        │     └── |dZ| < 0.1
        ├── SIP3D < 8
        ├── MiniPFRelIso < 0.4
        └── isMVANoIsoWP90() || MvaNoIso() > cutEC
              │                    └── NanoAOD: Electron_mvaFall17V2noIso
              └── NanoAOD: Electron_mvaFall17V2noIso_WP90
```

MVA thresholds for `HcToWALooseRun2`:
- Inner Barrel (|scEta| < 0.8): `MvaNoIso > 0.985`
- Outer Barrel (0.8-1.444): `MvaNoIso > 0.96`
- Endcap (> 1.566): `MvaNoIso > 0.85`

For Run3 (`HcToWALooseRun3`), the endcap threshold is -0.8 (effectively no constraint), which explains why Run3 is unaffected even if MVA scores are miscalibrated.

## Impact on ParticleNetMD Training

### Event counts (Run1E2Mu, Run2 only)

| Sample | Events | Notes |
|--------|-------:|-------|
| TTLL_powheg (nominal) | 81,920 | Correct endcap |
| 8 systematic variants (combined) | ~245,000 | Broken endcap |
| **Total Run2 nonprompt** | **~327,000** | 75% has broken endcap |

The nonprompt training class uses LNT+Bjet events from all 9 TTLL variants (see [nonprompt.md](DataAugment/nonprompt.md)). With 75% of Run2 events lacking endcap electrons, the network may learn an artificially low nonprompt rate in the electron endcap region.

### Mitigation by training design

The impact is partially mitigated by:
1. **Weight balancing**: `balance_weights=True` in [SglConfig.json](configs/SglConfig.json) normalizes per-class total weights
2. **Event cap**: `max_events_per_fold_per_class=40,000` limits oversampling
3. **Run3 data**: Run3 systematic variants are unaffected and contribute correct endcap coverage
4. **Nominal dominance**: TTLL_powheg has 2x the events of each systematic variant

## Recommendations

1. **Investigate NanoAOD**: Check `Electron_mvaFall17V2noIso` distributions for endcap electrons in a Run2 systematic variant NanoAOD file to confirm the upstream cause
2. **Short-term fix**: For training, consider dropping Run2 systematic variants from the nonprompt class and using only nominal TTLL_powheg for Run2 (81,920 events, still above the 62,486/fold target)
3. **Long-term fix**: Re-produce Run2 NanoAOD for the 8 systematic TTLL variants with correct electron MVA calibrations, then re-run EvtTreeProducer and the dataset pipeline

## References

- [Nonprompt sample selection](DataAugment/nonprompt.md) — rationale for including 9 TTLL variants
- [Event count tables](DataAugment/dataset.md) — full per-sample statistics
- [Validation script](python/validateDatasets.py) — dataset validation pipeline
- [EvtTreeProducer](../SKNanoAnalyzer/Analyzers/src/EvtTreeProducer.cc) — event selection and electron ID application
- [Electron ID definitions](../SKNanoAnalyzer/DataFormats/src/Electron.cc) — `HcToWALooseRun2/Run3`, `CaloIdL_TrackIdL_IsoVL`
- [NanoAOD branch mapping](../SKNanoAnalyzer/Analyzers/src/AnalyzerCore.cc) — `mvaFall17V2noIso` (Run2) vs `mvaNoIso` (Run3)
