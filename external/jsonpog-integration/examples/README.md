# Examples

Basic examples for standard usage of JSON tables with `correctionlib` are provided.

All examples are provided at least in Python, but some are also provided in C++. Both essentially do the same, with minor differences. The Python example can be directly executed in the command line, whereas the C++ example should first be compiled with `make`.

## BTV

Four use cases are shown
- fixedWP correction with mujets (here medium WP)
- fixedWP correction uncertainty (here tight WP and comb SF)
- shape correction SF 
- shape correction SF uncertainties

The results from C++ and Python examples differ in values, as the Python (C++) executable relies on NumPy / ROOT for the random number generation.

## Electrons

Examples for handling electron corrections are provided in Python. 
These examples demonstrate how to retrieve and apply scale factors (SFs) and corrections for different scenarios, including Run 2 UL and Run 3 campaigns.
The examples cover the following aspects.

### Electron ID Scale Factors

- Retrieve SFs for different sources and regimes (e.g., `RecoBelow20`, `RecoAbove20`; `Medium` for ID).
- Evaluate systematic variations.
- Example file: `electronExample.py`

### Electron HLT Scale Factors

- Retrieve efficiencies for HLT triggers with respect to a certain offline selection (e.g., `HLT_SF_Ele30_TightID`) in both data and MC.
- Retrieve the corresponding nominal SFs.
- Evaluate systematic variations.
- Example file: `electronHltExample.py`

### Electron Scale and Smearing Corrections

- Apply energy scale corrections for data and energy resolution (smearing) corrections for MC.
- Evaluate uncertainties for both scale and smearing corrections (to be applied to MC).
- Example file: `egmScaleAndSmearingExample.py`
- Please note: Usage of scale and smearing for photons is similar and we do not provide an explicit usage example.

## JERC

Six use cases are shown:
- Retrieve the jet energy correction at a single level (e.g. `L2Relative`) for AK4 and AK8.
- Retrieve the jet energy compound correction (all levels---this is probably what you need in most cases).
- Retrieve the jet energy uncertainty (e.g. `Total`---probably also what most analyses will need).
- Retrieve the jet energy resolution scale factor
- Retrieve the jet energy resolution 
- Retrieve the jet energy correction factor from the jet energy resolution smearing

## JMAR

Example file `jmarExample.py` contains examples of retrieving SFs in six scenarios:
- DeepAK8/ParticleNet tagging
- cut-based top tagging
- cut-based W tagging
- soft drop mass correction
- PU JetID
- Quark-Gluon tagging

Example file `jetidExample.py` contains examples of obtaining JetID decision flag.

Please note:\
JetID Jsons have been updated and ```chMultiplicity, neMultiplicity, multiplicity``` are now integer variables, reflecting the types in NanoAOD.\
However, these jsons will only work with ```correctionlib 2.6.0``` and later.\
If you receive an error similar to ```RuntimeError: std::get: wrong index for variant```, you are most likely using older version of ```correctionlib``` so please update to ```correctionlib 2.6.0+```.

## MET $\phi$

## Muons

Retrieving various SFs for Run 2 UL and Run 3 campaigns.

## Photons

Examples for handling photon corrections are provided in Python and C++. 
These examples demonstrate how to retrieve and apply scale factors (SFs) for different scenarios, including Run 2 UL and Run 3 campaigns.
The examples currently cover ID, conversion-safe electron veto (CSEV) and pixel veto (PixVeto) SFs.

### Photon SFs

- Retrieve SFs for different sources.
- Evaluate systematic variations.
- Example files: `photonExample.py`, `photonExample.C`.

## Tauons

