# Development Guide

## Contents

- [Repository Setup](#repository-setup)
- [Coding Standards](#coding-standards)
- [Pull Requests](#pull-requests)
- [GitHub Actions CI](#github-actions-ci)

---

## Repository Setup

Always start from a clean state:

```bash
git clone --recurse-submodules git@github.com:$GITACCOUNT/SKNanoAnalyzer.git
cd SKNanoAnalyzer
git remote add upstream git@github.com:CMSSNU/SKNanoAnalyzer.git

# Sync with upstream
git fetch upstream && git rebase upstream/main

# Create development branch
git checkout -b $DEVBRANCH
```

---

## Coding Standards

### Use Enum Classes (Not Magic Numbers/Strings)

Avoid string keys—they cause confusion (`"POG_TIGHT"` vs `"TIGHT"` vs `"T"`).

```cpp
// Good: enum class
enum class ElectronID {
    NOCUT,
    POG_VETO,
    POG_LOOSE,
    POG_MEDIUM,
    POG_TIGHT,
    POG_HEEP,
    POG_MVAISO_WP80,
    POG_MVAISO_WP90,
    POG_MVANOISO_WP80,
    POG_MVANOISO_WP90,
};

if (electron.PassID(ElectronID::POG_TIGHT)) { ... }
```

### Use Physical Constants

Use `PhysicalConstants.h` (included via `AnalyzerCore.h`) instead of hardcoding values:

```cpp
// Good
if (fabs(mass - Z_MASS) < 15) return false;

// Bad
if (fabs(mass - 91.1876) < 15) return false;
```

### Systematic Variations in Corrections

When adding corrections with systematic variations, use this signature pattern:

```cpp
// In MyCorrection.h
float GetYourCorrection(/* parameters */, const variation syst, const TString &source = "total");
```

The `variation` enum is defined in `MyCorrection.h`:

```cpp
enum class variation { nom, up, down };
```

Use `source = "total"` for combined uncertainties. See `GetBTaggingSF` for a reference implementation.

---

## Pull Requests

### Testing PRs from Others

```bash
# From your development branch
git fetch upstream pull/$PRNUMBER/head:pr$PRNUMBER
git checkout pr$PRNUMBER

# Example: test PR #16
git fetch upstream pull/16/head:pr16
git checkout pr16
```

### Opening a Pull Request

1. Push your branch to your fork
2. Open PR against `CMSSNU/SKNanoAnalyzer:main`
3. Ensure CI passes (see below)
4. Request review from maintainers

---

## GitHub Actions CI

Changes to CMSSNU/SKNanoAnalyzer trigger automatic builds. The CI requires SSH access to gitlab.cern.ch for the RoccoR submodule.

### Setting Up SSH Key for CI

1. Generate a dedicated key:

   ```bash
   ssh-keygen -t ed25519 -C "gitlab-ci-key" -f ~/.ssh/id_ed25519_gitlab_gha
   ```

2. Add **public key** to GitLab: [gitlab.cern.ch](https://gitlab.cern.ch) → Preferences → SSH Keys

3. Add **private key** to GitHub: Repository → Settings → Secrets and variables → Actions → New repository secret → Name: `SSH_PRIVATE_KEY`

> Note: `pull_request_target` triggers workflows from the upstream repo, so forked repositories don't need additional configuration.
