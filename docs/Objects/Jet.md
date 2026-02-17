# Jet

The `Jet` class represents reconstructed jets. Inherits from `Particle` → `TLorentzVector`.

## Contents

- [Class Hierarchy](#class-hierarchy)
- [Getting Jets](#getting-jets)
- [Identification](#identification)
- [B-Tagging](#b-tagging)
- [Energy Corrections](#energy-corrections)
- [Properties Reference](#properties-reference)
- [Run 2 vs Run 3](#run-2-vs-run-3)

---

## Class Hierarchy

```text
TLorentzVector
└── Particle        (charge, 4-momentum)
    └── Jet
```

**Files:** `DataFormats/include/Jet.h`, `DataFormats/src/Jet.cc`

---

## Getting Jets

### From AnalyzerCore

```cpp
// All jets
RVec<Jet> allJets = GetAllJets();

// With ID, pT, eta cuts
RVec<Jet> jets = GetJets("tight", 30., 2.4);

// Select from existing collection (enum-based)
RVec<Jet> selected = SelectJets(jets, Jet::JetID::TIGHT, 30., 2.4);

// Remove jets overlapping with leptons
RVec<Jet> cleanJets = JetsVetoLeptonInside(jets, electrons, muons, 0.4);
```

---

## Identification

### Enum-Based IDs (Recommended)

```cpp
enum class JetID {
    NOCUT,
    LOOSE,
    TIGHT,
    TIGHTLEPVETO,
    PUID_LOOSE,
    PUID_MEDIUM,
    PUID_TIGHT,
};

if (jet.PassID(Jet::JetID::TIGHT)) { ... }
```

### String-Based IDs

| String            | Description                       |
| ----------------- | --------------------------------- |
| `"loose"`         | Basic quality cuts                |
| `"tight"`         | Standard analysis quality         |
| `"tightLepVeto"`  | Tight + lepton veto               |
| `"puidLoose"`     | Pileup ID Loose (low-pT jets)     |
| `"puidMedium"`    | Pileup ID Medium                  |
| `"puidTight"`     | Pileup ID Tight                   |

### Direct ID Checks

```cpp
bool passTight    = jet.Pass_tightJetID();
bool passLepVeto  = jet.Pass_tightLepVetoJetID();
```

---

## B-Tagging

### Available Taggers

```cpp
namespace JetTagging {
    enum class JetFlavTagger {
        DeepJet,      // Legacy (Run 2/early Run 3)
        ParticleNet,  // Recommended for Run 3
        RobustParT,   // Alternative for systematics
    };
}
```

### Getting Tagger Results

```cpp
// B-tagging score
float bScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::ParticleNet);

// C-tagging scores (CvB, CvL)
auto [cvb, cvl] = jet.GetCTaggerResult(JetTagging::JetFlavTagger::ParticleNet);

// QvG score
float qvg = jet.GetQvGTaggerResult(JetTagging::JetFlavTagger::ParticleNet);
```

### Working Points

```cpp
// Get WP threshold
float btagWP = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParticleNet,
                                      JetTagging::JetFlavTaggerWP::Medium);

// Check if b-tagged
bool isBTagged = jet.GetBTaggerResult(JetTagging::JetFlavTagger::ParticleNet) > btagWP;
```

### Scale Factors

```cpp
// B-tagging SF
float btag_sf = myCorr->GetBTaggingSF(jet.Eta(), jet.Pt(), jet.hadronFlavour(),
                                       JetTagging::JetFlavTagger::ParticleNet,
                                       JetTagging::JetFlavTaggerWP::Medium,
                                       MyCorrection::variation::nom);

// B-tagging efficiency (from MC)
float btag_eff = myCorr->GetBTaggingEff(jet.Eta(), jet.Pt(), jet.hadronFlavour(),
                                         JetTagging::JetFlavTagger::ParticleNet,
                                         JetTagging::JetFlavTaggerWP::Medium,
                                         MyCorrection::variation::nom);
```

---

## Energy Corrections

### JES/JER Methods

```cpp
// Jet Energy Resolution
float jer = myCorr->GetJER(jet.Eta(), jet.Pt(), rho);

// JER Scale Factor
float jerSF = myCorr->GetJERSF(jet.Eta(), jet.Pt(), MyCorrection::variation::nom);

// JES correction
float jesSF = myCorr->GetJESSF(jet.GetArea(), jet.Eta(), jet.Pt(),
                                jet.Phi(), rho, runNumber);

// JES uncertainty
float jesUnc = myCorr->GetJESUncertainty(jet.Eta(), jet.Pt(),
                                          MyCorrection::variation::up, "Total");
```

### Correction Properties

```cpp
float rawPt      = jet.GetRawPt();           // Uncorrected pT
float originalPt = jet.GetOriginalPt();      // Original stored pT
float rawFactor  = jet.j_rawFactor;          // Raw factor

// Regression corrections
float pnetCorr = jet.j_PNetRegPtRawCorr;     // ParticleNet regression
float bRegCorr = jet.j_bRegCorr;             // b-jet regression
float cRegCorr = jet.j_cRegCorr;             // c-jet regression
```

---

## Properties Reference

### Energy Fractions

```cpp
float chHEF  = jet.chHEF();      // Charged hadron
float neHEF  = jet.neHEF();      // Neutral hadron
float chEmEF = jet.chEmEF();     // Charged EM
float neEmEF = jet.neEmEF();     // Neutral EM
float muEF   = jet.muEF();       // Muon
```

### Multiplicities

```cpp
short nConst  = jet.nConstituents();    // Total constituents
short nElec   = jet.nElectrons();       // Electrons
short nMuons  = jet.nMuons();           // Muons
short nSVs    = jet.nSVs();             // Secondary vertices

// NanoAODv13
uchar chMult  = jet.chMultiplicity();   // Charged hadron multiplicity
uchar neMult  = jet.neMultiplicity();   // Neutral hadron multiplicity
```

### Flavor Information (MC)

```cpp
int partonFlav = jet.partonFlavour();   // Parton-based flavor
int hadronFlav = jet.hadronFlavour();   // Hadron-based flavor
int genJetIdx  = jet.genJetIdx();       // Matched GenJet index
```

### Geometry

```cpp
float area = jet.GetArea();             // Jet area (for JEC)
float mass = jet.GetM();                // Jet mass
```

---

## Run 2 vs Run 3

| Aspect         | Run 2 (NanoAODv9)        | Run 3 (NanoAODv13)             |
| -------------- | ------------------------ | ------------------------------ |
| PU Algorithm   | CHS                      | PUPPI                          |
| JEC Levels     | L1+L2+L3+residual        | L2+L3+residual                 |
| Jet ID         | NanoAOD bits             | Custom implementation          |
| B-Tagger       | DeepJet primary          | ParticleNet / RobustParT       |
| Eras           | 2016pre/post, 2017, 2018 | 2022, 2022EE, 2023, 2023BPix   |

> **Note:** Run 3 uses custom JetID due to NanoAODv12 bug.

---

## References

- [CMS Jet ID 13.6 TeV](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV)
- [JEC/JER Recommendations](https://cms-jerc.web.cern.ch/Recommendations/)
- [Run 3 Jet Performance](http://cds.cern.ch/record/2911750)
- [BTV POG Recommendations](https://btv-wiki.docs.cern.ch/)
