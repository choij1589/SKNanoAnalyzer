# SKNanoAnalyzer
---

## THIS IS DEVELOPMENT VERSION!!

## Introduction
- If you want to engage in the development and commit to the repository, please read [Development Guide](docs/DevelopmentGuide.md) first.
- For setting up the environment and starting the analysis, please refer [Getting Started](docs/GettingStarted.md).

## To do
### Assigned
- [x] Add new features to ExampleRun.cc and make them manual: **Yeonjoon**
- [ ] Check the update in jsonpog-correction when sourcing setup.sh: **Chihwan**
- [x] Clean up gv0: **Jin**
- [ ] Compare the sample list in SKFlat with the necessary samples, and bring them if necessary, check xsec, etc.: **Taehee+Eunsu+Jin**
- [x] Update the CheckStorages.sh to respect the folder structure under gv0: **Jin**
- [ ] FatJet and Tau: **Youngwan**
- [ ] SKNano.py skimming test **Yeonjoon**
- [ ] Auto-Validation: Check the consistency of the output for a specific NanoAOD file: **Jin**
- [x] Check for CVMFS environment compatability **deprecated**
- [x] Update MeasureTaggingEff.cc **Yeonjoon**
- [ ] Rochester Correction Follow-up **Ask Hyonsan or Won(Contact: Chihwan)**
- [ ] Update the network to 10G? (later)

For DY, check the before/after of correction
- [ ] Leptons **Jihoon**
- [ ] Jets(B-tagging) **Taehee**

## Make Tutorial
- [ ] From make Training tree, Export ONNX, Use this model in the Analyzer : **Yeonjoon(Review: Youngwan)**
- [ ] Simple DY tutorial : **Eunsu(Review: Jin)**
- [ ] tt reconstruction using Kinematic Fit : **Minseok(Review: Byunghun)**
- [ ] tt reconstruction using KinFit Library : **Yeonjoon**
- [ ] NLO vs LO comparison in MadGraph Tutorial: **Jin**
- [ ] Tau+Fatjet Tutorial : **Youngwan**

- [ ] At least try to install until next meeting : ***All***


## Useful links
- [MiniAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD)
- [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD)
- [NanoAOD-Tools (for postproc)](https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAODTools)
- [SNU_SKFlatProduction](https://docs.google.com/spreadsheets/d/12SfHxUR_1_-13Tzt6mI0slnenQeeCT40zDdggx_KPR8/edit?pli=1#gid=1739893496)
- [NanoAOD wiki@gitlab](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/home)
- [Run2 Summer20 UL](https://docs.google.com/spreadsheets/d/1ABl2p2uwr2EfEbolBEVNcKb_fIXigYY9sqCRT8XIi1Q/edit#gid=1318927481)
- [Run3 Common Bkg Sample Preparation](https://docs.google.com/spreadsheets/d/1xEbHtxzJVaWJpxDrJ7tSpxcVQMRIFbZeoeO7davNFMA/edit#gid=1396344675)
- [SKFlatMaker](https://github.com/CMSSNU/SKFlatMaker)
- [SKFlatAnalyzer](https://github.com/CMSSNU/SKFlatAnalyzer)
- [SKFlat Manual](https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf)
- [NanoAOD Branches](https://cms-nanoaod-integration.web.cern.ch/autoDoc/)
- [PdmV Run3 DATA full table](https://pdmv-pages.web.cern.ch/run_3_data/full_table.html)
- [CMS 13.6TeV SM cross sections](https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV)
- [LHC ttbar NNLO+NNLL cross sections](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO)
- [LHC Single Top NNLO+NNLL cross sections](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec)
- [LHC Higgs Cross Section](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/HiggsXSBR)
- [Rules for Run3 2024 dataset names](https://cms-pdmv.gitbook.io/project/mccontact/rules-for-run3-2024-dataset-names)
- [Rules for Run3 dataset names](https://cms-pdmv.gitbook.io/project/mccontact/rules-for-run3-dataset-names)
