# Jet Veto Map Implementation Notes

## Overview

This document describes the jet veto map implementation in the DiLepton analysis and a known issue with the Run2 efficiency calculation.

## Run2 vs Run3 Strategy

### Run2: Jet-level veto
- Individual jets in problematic detector regions are removed
- Other jets in the event are kept
- Efficiency metric: `(jets after veto) / (jets before veto)`

### Run3: Event-level veto
- If ANY jet falls in the veto map region, the entire event is rejected
- Efficiency metric: `(events after veto) / (events before veto)`
- Calculated from `ALL/Central/cutflow` bins 2 (NoiseFilter) → 3 (VetoMap)

## Known Issue: Run2 Efficiency Calculation

### Problem

The current Run2 jet veto efficiency appears artificially low (~60% for data, ~80% for MC) and shows **uniform** rejection across all eta-phi regions, which is inconsistent with a localized veto map.

### Root Cause

In `AnalyzerCore::PassVetoMap(const Jet&, ...)` (line 1162-1176 of AnalyzerCore.cc):

```cpp
bool pass_loose_selection = jet.Pt() > 15.;
pass_loose_selection = pass_loose_selection && jet.PassID(Jet::JetID::TIGHT);
pass_loose_selection = pass_loose_selection && (jet.Pt() > 50. || jet.PassID(Jet::JetID::PUID_LOOSE));
for (const auto &muon: AllMuons){
    pass_loose_selection = pass_loose_selection && (jet.DeltaR(muon) > 0.2);
}
bool pass_veto_map = pass_loose_selection && (!myCorr->IsJetVetoZone(...));
return pass_veto_map;
```

The function returns `false` (jet removed) if `pass_loose_selection` is false, **even if the jet is NOT in a veto zone**.

The muon overlap check (`DeltaR(muon) > 0.2`) causes uniform rejection across all eta-phi regions.

### Current Histogram Filling in DiLepton.cc

```cpp
// BeforeJetVeto: filled with all tightJets (includes jets near muons)
FillJetEtaPhi2D(tightJets, 1.0, "BeforeJetVeto");

// PassVetoMap removes jets near muons AND jets in veto zones
for (const auto &jet: tightJets)
    if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.emplace_back(jet);

// AfterJetVeto: excludes jets near muons
FillJetEtaPhi2D(tightJets, 1.0, "AfterJetVeto");
```

The efficiency mixes two effects:
1. **Muon overlap removal** (~30-40%, uniform across eta-phi)
2. **Actual veto map removal** (small, localized to specific regions)

### Why This is Acceptable (for now)

Jets overlapping with muons would be removed later anyway by `JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4)`. So from the final analysis perspective, the result is the same.

However, for a **pure measurement** of the jet veto map efficiency, this is misleading.

### Suggested Fix for Future

To measure true jet veto map efficiency, change the order in `DiLepton.cc`:

```cpp
// 1. Apply muon overlap first (same criteria as PassVetoMap loose selection)
RVec<Jet> tightJets_baseline;
for (const auto &jet: tightJets) {
    bool pass_loose = (jet.chEmEF() + jet.neEmEF() > 0.9);  // EM jets bypass
    if (!pass_loose) {
        pass_loose = true;
        for (const auto &muon: allMuons)
            pass_loose = pass_loose && (jet.DeltaR(muon) > 0.2);
    }
    if (pass_loose) tightJets_baseline.emplace_back(jet);
}

// 2. Fill BEFORE with baseline (muon overlap already applied)
FillJetEtaPhi2D(tightJets_baseline, 1.0, "BeforeJetVeto");

// 3. Apply ONLY veto map check
RVec<Jet> tightJets_vetoMap;
for (const auto &jet: tightJets_baseline)
    if (!myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap"))
        tightJets_vetoMap.emplace_back(jet);

// 4. Fill AFTER
FillJetEtaPhi2D(tightJets_vetoMap, 1.0, "AfterJetVeto");
```

Alternatively, fix `PassVetoMap` to return `true` for jets failing loose selection:

```cpp
if (!pass_loose_selection) return true;  // Don't apply veto map, but keep the jet
return !myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory);
```

## Run3 Implementation: Verified Correct

The Run3 event-level veto in `AnalyzerCore::PassVetoMap(const RVec<Jet>&, ...)` is correctly implemented:

```cpp
bool AnalyzerCore::PassVetoMap(const RVec<Jet> &AllJets, const RVec<Muon> &AllMuons, ...) {
    if (! (Run == 3)) return true;  // Only for Run3

    // Select jets meeting loose criteria
    RVec<Jet> this_jet = SelectJets(AllJets, Jet::JetID::TIGHT, 15., 5.0);
    this_jet = JetsVetoLeptonInside(this_jet, empty_electrons, AllMuons, 0.2);
    for(const auto &jet: this_jet){
        if(jet.chEmEF() + jet.neEmEF() < 0.9) selected_jets.push_back(jet);
    }

    // If ANY jet is in veto zone, reject entire event
    for(const auto &jet: selected_jets){
        if(myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory)) return false;
    }
    return true;
}
```

## Plotting Script

The plotting script `python/plotJetVetoMap.py` generates:
- 2D eta-phi distributions (before/after for Run2, passed for Run3)
- Efficiency summary (printed and saved to `efficiency.json`)

Usage:
```bash
python python/plotJetVetoMap.py --era 2018 --channel DIMU
python python/plotJetVetoMap.py --era 2022 --channel DIMU
```

Output location: `plots/{era}/{channel}/JetVetoMap/`
