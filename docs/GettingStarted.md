# Getting Started

## Important - Skim First

The main bottleneck on SNU Server is I/O. Most events fail baseline selection (~80-90%), causing unnecessary I/O without skimming. **Always skim samples before running full analysis.**

After skimming, submit all matching samples with one command:

```bash
SKNano.py -a MyAnalyzer -i 'Skim_DiLepton_*' -e 2022 -n 10
```

## Contents

- [Environment Setup](#environment-setup)
- [Installation](#installation)
- [Job Submission](#job-submission)
- [Making Sample Lists](#making-sample-lists)
- [Skimming Mode](#skimming-mode)
- [Telegram Notifications](#telegram-notifications)
- [Appendix: macOS Setup](#appendix-macos-setup)

---

## Environment Setup

### Package Manager (micromamba recommended)

Micromamba is faster than conda for environment solving:

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

**Important:** Install packages to `/data6/Users/$USER/`, not home directory (inaccessible from worker nodes).

### Create Environment

#### Option 1 - From template (recommended)

```bash
cd $SKNANO_HOME/docs
# Edit Nano.yml: change prefix to your path
micromamba env create -f Nano.yml
```

#### Option 2 - Manual setup

```bash
micromamba create -n nano python=3.12 root=6.34.04 -c conda-forge
micromamba activate nano

# Required packages
micromamba install onnxruntime-cpp correctionlib boost-cpp -c conda-forge

# Optional: ML packages
pip install torch==2.4.1 --index-url https://download.pytorch.org/whl/cu121
pip install torch_geometric numpy pandas matplotlib scipy scikit-learn
```

### SSH Key for CERN GitLab

Required for RoccoR submodule:

```bash
ssh-keygen -t ed25519 -C "your.email@cern.ch"
```

Add public key at: [gitlab.cern.ch](https://gitlab.cern.ch) → Preferences → SSH Keys

---

## Installation

### Clone Repository

```bash
git clone --recurse-submodules git@github.com:$GITACCOUNT/SKNanoAnalyzer.git
cd SKNanoAnalyzer
git remote add upstream git@github.com:CMSSNU/SKNanoAnalyzer.git
```

### Configure

```bash
cp config/config.default config/config.$USER
```

Edit `config/config.$USER`:

| Variable            | Description                                                                                    |
| ------------------- | ---------------------------------------------------------------------------------------------- |
| `SYSTEM`            | `redhat` or `osx`                                                                              |
| `PACKAGE`           | `conda` or `mamba`                                                                             |
| `SKNANO_HOME`       | Path to repository (Optional)                                                                  |
| `SKNANO_RUNLOG`     | Path for run logs (Optional)                                                                   |
| `SKNANO_OUTPUT`     | Path for output files (Optional)                                                               |
| `TOKEN_TELEGRAMBOT` | Telegram bot token (optional)                                                                  |
| `USER_CHATID`       | Your Telegram chat ID (optional)                                                               |
| `SINGULARITY_IMAGE` | Container image path (optional, e.g., `/data9/Users/choij/Singularity/images/private-el9.sif`) |

### Build

```bash
source setup.sh        # Required every session (auto-installs LHAPDF/libtorch on first run)
./scripts/build.sh     # Full build
```

### Singularity Support

For batch jobs using Singularity, build inside the container:

```bash
singularity exec $SINGULARITY_IMAGE bash -c "source setup.sh && ./scripts/build.sh"
```

### Verify Installation

```cpp
// ROOT
root -l
Particle *p = new Particle;
p->SetPtEtaPhiM(30, 2.1, 1.3, 0.1);
p->Print()
```

```python
# Python
from ROOT import Particle
p = Particle()
p.SetPtEtaPhiM(30, 2.1, 1.3, 0.1)
p.Print()
```

---

## Job Submission

Submit jobs to HTCondor using `SKNano.py`:

```bash
SKNano.py -a AnalyzerName -i SampleName -e era -n files_per_job
```

### Options

| Flag              | Description                                                                        |
| ----------------- | ---------------------------------------------------------------------------------- |
| `-a`              | Analyzer name                                                                      |
| `-i`              | Sample pattern (regex supported, e.g., `'DYJets*'`)                                |
| `-e`              | Era: `2022`, `2022EE`, `2023`, `2023BPix`, `2016preVFP`, etc. (comma-separated)    |
| `-n`              | Files per job (positive) or total jobs (negative, e.g., `-n -100` = 100 total)     |
| `-r`              | Run period: `Run2` or `Run3` (overrides `-e`)                                      |
| `--reduction N`   | Process 1/N events                                                                 |
| `--memory`        | Job memory (default: 2GB)                                                          |
| `--ncpu`          | CPUs per job (default: 1)                                                          |
| `--userflags`     | Custom flags (comma-separated)                                                     |
| `--batchname`     | Custom batch name                                                                  |
| `--skimming_mode` | Enable skimming output                                                             |
| `--no_exec`       | Generate DAG without submitting                                                    |

### Examples

```bash
# All samples starting with "ST"
SKNano.py -a Vcb_FH -i 'ST*' -e 2022EE -n 10

# Multiple eras
SKNano.py -a DiLepton -i DYJets -e 2022,2022EE -n 50

# With reduction for testing
SKNano.py -a DiLepton -i TTLJ_powheg -e 2022 -n 10 --reduction 100
```

---

## Making Sample Lists

After saving NanoAOD samples to `/gv0`:

1. **Update sample metadata** in `data/$ERA/Sample/CommonSampleInfo.json`

2. **Generate path info:**

   ```bash
   ./scripts/MakeSamplePathInfo.py --era $ERA
   ```

3. **Calculate effective luminosity:**

   ```bash
   SKNano.py -a GetEffLumi -i $SAMPLENAME -e $ERA -n 10
   ```

4. **Update metadata with results:**

   ```bash
   ./scripts/parseEffLumi.py --era $ERA
   ```

---

## Skimming Mode

Skimming creates reduced samples at `$SKNANO_RUN[2,3]_NANOAODPATH/Era/[Data,MC]/Skim/$USERNAME/`.

### Submit Skimming Job

```bash
# Explicit flag
SKNano.py -a AnalyzerName -i DYJets -e 2022 -n 10 --skimming_mode

# Or use Skim_ prefix (prompts for confirmation)
SKNano.py -a Skim_DiLepton -i DYJets -e 2022 -n 10
```

### Use Skimmed Samples

After skimming completes:

```bash
SKNano.py -a MyAnalyzer -i Skim_DiLepton_DYJets -e 2022 -n 10
```

**Warning:** Do not submit multiple skimming DAGs simultaneously—they modify `skimTreeInfo.json` sequentially.

---

## Telegram Notifications

1. **Create bot:** Message `@BotFather` on Telegram, send `/newbot`, follow instructions

2. **Get token:** Save the token to `TOKEN_TELEGRAMBOT` in your config

3. **Get chat ID:** Send any message to your bot, then visit:

   ```text
   https://api.telegram.org/bot[YOUR_TOKEN]/getUpdates
   ```

   Find `"from":{"id":YOUR_CHAT_ID...` in the response

4. **Save chat ID** to `USER_CHATID` in your config

---

## Appendix: macOS Setup

macOS works for development but has GPU/ML limitations.

```bash
# Install micromamba
brew install micromamba

# Create environment
micromamba create -n nano python=3.12 root=6.34.04 -c conda-forge
micromamba activate nano
micromamba install correctionlib onnxruntime-cpp boost-cpp -c conda-forge

# ML packages
pip install torch==2.4.1 torch_geometric
pip install --no-build-isolation git+https://github.com/pyg-team/pyg-lib.git
pip install --no-build-isolation torch_scatter torch_sparse torch_cluster torch_spline_conv
pip install numpy pandas matplotlib scipy scikit-learn networkx
```

### Building ROOT from Source (if conda version fails)

```bash
micromamba activate nano
cd ~/Downloads
git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git root_src

mkdir -p ~/mamba/root_build ~/mamba/root_install
cd ~/mamba/root_build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/mamba/root_install -Dbuiltin_glew=ON $HOME/Downloads/root_src
cmake --build . --target install -j8

# Cleanup
rm -rf ~/Downloads/root_src ~/mamba/root_build

# Link onnxruntime
ln -s $HOME/mamba/envs/nano/lib/libonnxruntime.*.dylib $HOME/mamba/root_install/lib/

# Add to shell
echo "source $HOME/mamba/root_install/bin/thisroot.sh" >> ~/.zshrc
```

Tested on: M4 Mac Mini, macOS Tahoe 15.0, Python 3.12, ROOT 6.34.04
