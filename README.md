# SKNanoAnalyzer
---

## THIS IS DEVELOPMENT VERSION!!

## Project Updates
[HedgeDoc](https://demo.hedgedoc.org/VrWRIlceTjO9SPOVLKUNVA?view)

## Notes
- ROOT v6.32.xx have memory leak issue while reading root files with <= v6.30.xx -> Avoid using ROOT >= v6.32.xx temporariliy

## To do
- Make ExampleRun runable -> done
- DataFormats -> Electron, Muon, LHE(Jin), Jet/GenJet (Yeonjun), Tau/FatJet(Youngwan), Gen(Taehee)
- TODO: start validation task, check installation of correctionlib, private NanoAOD generation
- external dependencies: ~~LHAPDF~~ / GEScaleSyst(?) / ~~CorrectionLib~~

## How to
### Preliminary Setup
Here is an example to setup the environment using conda.
```bash
# create conda environment
conda create -n nano python=3.11 root=6.30.04 -c conda-forge
conda activate nano

# Install onnxruntime-cpp
conda install onnxruntime-cpp

# Install correctionlib
pip install correctionlib
```

### Installation
Recommend to fork the repo to your account.
```bash
# Clone the repository
# jsonpog-integration should be cloned as a submodule if you don't use cvmfs
git clone --recurse-submodules git@github.com:$GITACCOUNT/SKNanoAnalyzer.git

# Add to your remote repo
git remote add upstream git@github.com:CMSSNU/SKNanoAnalyzer.git
git checkout $DEVBRANCH

# create config file
cp config/config.default config/config.$USER
# edit the configuration!

# first time setup
source setup.sh    # you have to do this every new session
```

### Testing PRs
Here is a recommended way to test PRs.
```bash
# Always start from the clean state
git clone --recurse-submodules git@github.com:$GITACCOUNT/SKNanoAnalyzer.git
cd SKNanoAnalyzer
git remote add upstream git@github.com:CMSSNU/SKNanoAnalyzer.git

# fetch the PR
git fetch upstream pull/$PRNUMBER/head:pr$PRNUMBER
git checkout pr$PRNUMBER

# for example, fetching PR 16
git fetch upstream pull/16/head:pr16
git checkout pr16
```

### Compilation
We are using cmake for the default compiling management.
Use scripts/build.sh for clean compilation.
```bash
./scripts/build.sh
```

Or you can do it manually
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$SKNANO_HOME ..
make
make install
cd $SKNANO_HOME
```

## Check modules
Every module(or class) can be imported both in ROOT and python
```cpp
root -l
Particle *p = new Particle;
p->SetPtEtaPhiM(30, 2.1, 1.3, 0.1);
p->SetCharge(1);
p->Print()
```

```python
python
from ROOT import Particle
p = Particle()
p.SetPtEtaPhiM(30, 2.1, 1.3, 0.1)
p.SetCharge(1)
p.Print()
```

For testing other modules and analyzers, check scripts/test.py

## About LHAPDFs
For using LHAPDFHandler and PDFReweight classes, two possible options
1. install lhapdf manually.
```bash
./scripts/install_lhapdf.sh
```
It would be run automatically for the first time setup.
2. use lhapdf from cvmfs

## About correctionlibs
In the config/config.$USER file, there is an option to choose bewteen conda and cvmfs. When configuring your environment with conda, at least ROOT and correctionlibs should be installed:
```bash
# example
conda env create -n nano python=3.11 root=6.30.04 -c conda-forge
conda activate nano
pip install correctionlib
```
## How to Submit the job

Jobs can be submitted to htcondor using SKFlat.py
```bash
SKFlat.py -a AnalyzerName -i SamplePD -n number of jobs -e era
```

Basic usage is as aboves. There are some additional options for the submission:
- -i: You can pass the sample PD using this option. This option supports the basic regex, thus,
  ```bash
  SKFlat.py -a Vcb_FH -i 'ST*' -n 100 -e 2022EE
    ```
    will submit the jobs for all the samples starting with 'ST' in the sample list. Thus, please be careful when you adding the new samples to the sample list.
- -n: Number of jobs to submit. If you want to submit 100 jobs, you can use -n 100.
  You can also choose to set the ***number of files for each job***. To do this, pass this argument as negative value. For example, -n -10 will submit 10 files per job. For example, if *TTLJ_powheg* sample has 1700 files, 
  ```bash
    SKFlat.py -a Vcb_FH -i TTLJ_powheg -n -10 -e 2022EE
    ```
    will submit 170 jobs with 10 files each to cluster.
- -e: Era of the sample. You can also pass the multiple eras using comma. For example, -e 2022EE,2023 will submit the jobs for the samples in 2022EE and 2023 era.
- -r: This argument set the run. choose Run2 or Run3. This option overrids the -e option.
- --reduction: perform reduction by factor of input number. 
- --memory: set the memory for the job. Default is 2GB.
- --ncpu: set the number of cpus for the job. Default is 1.
- --userflags: set the user flags for the job. Default is empty. to set multiple flags, use comma. e.g. --userflags flag1,flag2
- --batchname: set the batch name for the job. Default is the analyzer name_userflags.
- --skimming_mode: by passing this flag, SKFlat.py will submit the jobs for the skimming mode. Detailed information as follows.

### Skimming mode
By passing --skimming_mode, SKFlat.py will submit the jobs for the skimming mode. In this mode, the jobs will create the output in `$SKNANO_RUN[2,3]_NANOAODPATH/Era/[Data,MC]/Skim/$USERNAME` directory, Instead of submit hadd layer in DAG, *PostProc* layer will add in the DAG. 
If your analyzer has name that starts with "Skim_", you will be asked to be enable the skimming mode. If you choose to enable the skimming mode, then skimming mode will be activated. Of course you can manually activate the skimming mode by passing --skimming_mode flag.
PostProc layer will create the Skimmed sample folder and will creat the skimTreeInfo.json file and dedicated json that saves the information of the skimmed samples under $SKNANO_DATA/era/Sample/Skim directory.
Each postproc job modify the skimTreeInfo.json sequentially, so ***DO NOT SUBMIT THE DAG CLUSTER THAT DO SKIMMING.*** After all the postproc jobs are done, you can submit the new jobs that using skimmed sample. A prefix of Skim_AnalyzerName will be added to the output file name.
```bash
SKFlat.py -a AnalyzerName -i DYJets -n -1 -e era --skimming_mode
```
or 
```bash
SKFlat.py -a Skim_AnalyzerName -i DYJets -n -1 -e era
```
Will create the Skim_AnalyzerName_DYJets (if you choose to answer "y" in latter one).
Then you can submit the jobs by
```bash
SKFlat.py -a AnalyzerName -i Skim_AnalyzerName_DYJets -n -1 -e era
```


## To do
- Make ExampleRun runable -> done
- DataFormats -> Electron, Muon, LHE(Jin), Jet/GenJet (Yeonjun), Tau/FatJet(Youngwan), Gen(Taehee)
- TODO: start validation task, check installation of correctionlib, private NanoAOD generation
- external dependencies: LHAPDF / GEScaleSyst(?) / CorrectionLib

## Note

# SKFlat - Run3 Update

## General
- Scope of the update: Run3 only
- Baseline: Central NanoAOD → postproc (Private Skim + add branch)
- Any missing variables in NanoAOD? Need processing from MiniAOD?
- Private NANO: MiniAOD $\to$ NanoAOD (JME Private) 
- Jihoon will follow (SKFlatMaker): need time

## Object updates
- Event(Trigger, MET): trigger list 조사 / json? 영완
- DataStream 이름 조사 - 영완
- Muon / Electron: 진 / 영완 (Correction updated)
- Tau / FatJet: 영완 믿습니다(brave student) (updated)
- JetTaggingParameter (updating): 영완
- Jet / GenJet: 연준 / 지훈 (updated)
- GEN: 태희 / 진 (updating)
- Photon: 지훈(neeeww) (updating)
    - isolation value만
SNUCMSSW? 지훈(neeeww)
- systematic variation update

## SKFlat.py refactoring
- Run3용으로 update
- Done: Dagman framwork, Telegram notification, No Data-period support, auto iterating submission

### Validation task
- Drell Yan Peak 보기
- Scale Variation
- LHAPDF reweight validation
- GEN + GENJET validation: flag가 맞는지 확인하는 정도


## To do
- [ ]  Make a list of Data / MC samples and check the replicas [google docs](https://docs.google.com/spreadsheets/d/1fyk_y6L2OoVgDvDmSCGPMKjQEKwaNSo0Ji2NXtNKRf8/edit?usp=sharing)
- [ ]  Postproc development → survey for objects / branches / corrections / scale factors -> after sample list finalized
- [ ]  Make Baseline 
- support for local run
- survey for /gv0 storage (not using skim files)
- 오브젝트 가져올때 Skim된 Nano에서 branch가 drop되어있는지 확인(e.g. flags) -> Youngwan will follow
- SNUCMS에 github repo 파기


RDataframe / uproot (memory issue) / SKFlat
SKFlat(base) -> std::vector<Muon>
             -> ROOT::RVec<Muon>
Correctionlib -> only python binding -> Youngwan will follow
    
- history fill -> hash map in C++?
    

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
