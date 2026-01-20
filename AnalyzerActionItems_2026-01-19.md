# Tri-Lepton Analyzer Action Items

**Date:** 2026-01-19
**Purpose:** Consolidated list of analyzer-level changes required before re-running

---

## Code Changes Required

### CRITICAL

- [x] **PileupJetID Scale Factors** (JME)
  - Apply PileupJetID efficiency scale factors to MC
  - Follow JME POG recommendations for SF application
  - Add PileupJetID uncertainty as systematic variation in templates

### HIGH

- [x] **HEM Electron Veto - 2018 Only** (Reza)
  - Extend HEM veto to reject events with electrons in HEM-affected region
  - Region: η < -1.25, -1.62 < φ < -0.82 (matches existing IsHEMElectron implementation)
  - Apply to both data and MC for 2018
  - Check the yield impact on 2018 data and MC
  - We will check the impact on the fake rate measurement for V4 -> V5 update
  - **IMPLEMENTED:** vetoHEM=true passed to SelectElectrons when DataYear==2018

- [x] **Cone-Corrected pT Fix** (Marina)
  - Fix cone-corrected pT implementation in nonprompt lepton estimation
  - Affects Run-3 normalisation in closure tests
  - Re-produce closure plots after fix

---

## New Histograms/Outputs

### HIGH

- [x] **2D Jet η-φ Distributions** (JME)
  - Di-Lepton Control Regions
  - Produce before jet veto map application
  - Produce after jet veto map application
  - For each era: 2016preVFP, 2016postVFP, 2017, 2018, 2022, 2022EE, 2023, 2023BPix
  - Purpose: Confirm HEM, EE+, BPix corrections are properly applied
  - **IMPLEMENTED:** FillJetEtaPhi2D helper added to AnalyzerCore with calls in DiLepton.cc
    - Run 2: BeforeJetVeto_Run2, AfterJetVeto_Run2 stages in defineObjects
    - Run 3: PassedEventVeto_Run3 stage in executeEvent

- [x] **Jet Veto Yield Loss Tracking** (Marina)
  - Tri-Lepton Signal Regions and TT+Z Control Regions
  - Track yield before and after jet veto maps
  - Output percentage loss for data and signal separately
  - Separate Run-2 (jet removal) vs Run-3 (event removal) implementations

### MEDIUM

- [x] **Low-LR Validation Region Plots** (Marina)
  - Produce data/MC comparison in region below ParticleNet LR cut
  - For on-Z mass points (80 < mA < 100 GeV)
  - Verify signal contamination is minimal

---

## Systematic Uncertainty Changes

### CRITICAL

- [x] **PileupJetID Uncertainty** (JME)
  - Add as nuisance parameter in fit
  - Propagate to signal and background templates

### HIGH

- [x] **JES Nuisance Parameter Naming** (JME)
  - Rename "JetEN" to standard CMS convention (e.g., CMS_scale_j)
  - Consider year-decorrelated scheme if post-fit shows significant constraint

- [x] **B-tagging Nuisance Verification** (BTV)
  - Verify 10 nuisance parameters correctly implemented for Run 2:
    - 8 uncorrelated: `btagSFbc_uncorr_YEAR`, `btagSFlight_uncorr_YEAR` (x4 years)
    - 2 correlated: `btagSFbc_corr`, `btagSFlight_corr`
  - Provide impact breakdown table

---

## Pre-Run Checklist

```
Before submitting jobs:
[x] PileupJetID SF code implemented
[x] PileupJetID systematic variation added
[x] HEM electron veto added for 2018
[x] Cone-corrected pT bug fixed
[x] 2D jet η-φ histograms added
[x] Jet veto yield counters added
[x] JES nuisance parameter renamed
[x] B-tagging nuisance parameters verified
```

---

## Post-Run Validation

```
After jobs complete:
[ ] Check 2D η-φ plots show correct veto regions
[ ] Verify jet veto yield loss is reasonable
[ ] Compare 2018 yields before/after electron HEM veto
[ ] Check closure plots with cone-corrected pT fix
[ ] Verify PileupJetID uncertainty appears in impact plots
[ ] Confirm JES naming in datacards
[ ] Confirm 10 b-tagging nuisances in datacards
```

---

## Reference: Review Sources

| Item | Reviewer | File |
|------|----------|------|
| PileupJetID | JME | ReviewerResponses_JME_BTV_2026-01-19.md |
| Jet veto documentation | JME | ReviewerResponses_JME_BTV_2026-01-19.md |
| JES clarification | JME | ReviewerResponses_JME_BTV_2026-01-19.md |
| B-tagging verification | BTV | ReviewerResponses_JME_BTV_2026-01-19.md |
| HEM electron veto | Reza | ReviewerResponses_Reza_2026-01-19.md |
| Cone-corrected pT | Marina | ReviewerResponses_Marina_2026-01-10.md |
| Jet veto yield loss | Marina | ReviewerResponses_Marina_2026-01-10.md |
| Low-LR validation | Marina | ReviewerResponses_Marina_2026-01-10.md |
