import numpy as np
import awkward as ak
import correctionlib
import gzip
import json
import os

###########################################
# Create example electron variables (toy) #
###########################################

# Create two large samples of one electron per event for both data and MC,
# where pT is sampled from a Gaussian distribution.
# This is not realistic but serves as an example, where we can check that scale and smearing work as intended.
# seedGain, supercluster eta and r9 are kept constant across all events for simplicity.
# Please note: The usage of photons is the same as for electrons, but the correction values are slightly different.
N = 1000000

rng = np.random.default_rng(seed=125)
data_pt = rng.normal(loc=50.0, scale=1.0, size=N) # Note np language: loc is the location (mean), scale is the standard deviation (sigma)
# In our language, scale affects the mean of the distribution, and smearing affects the standard deviation.
data_electrons = ak.Array({
    "seedGain": np.full(N, 12),
    "pt": data_pt,
    "ScEta": np.full(N, 0.0),
    "r9": np.full(N, 0.99)
})
data_run = np.full(N, 324219) # example run number for data

mc_pt = rng.normal(loc=50.0, scale=1.0, size=N) # In practice, the mean and sigma of the MC pT distribution would be different.
mc_electrons = ak.zip({
    "seedGain": np.full(N, 12),
    "pt": mc_pt,
    "ScEta": np.full(N, 0.0),
    "r9": np.full(N, 0.99)
})
mc_run = ak.Array(np.ones(N, dtype=int))  # run for each MC event (also set to 1 in nanoAOD)

###########################################
# Load the correctionlib JSON corrections #
###########################################

path_ETDependent = "../POG/EGM/2022_Summer22/electronSS_EtDependent.json"

if not os.path.exists(path_ETDependent):
    os.system(f"gunzip -k {path_ETDependent}.gz")

cset = correctionlib.CorrectionSet.from_file("../POG/EGM/2022_Summer22/electronSS_EtDependent.json")

scale_evaluator = cset.compound["EGMScale_Compound_Ele_2022preEE"]
smear_and_syst_evaluator = cset["EGMSmearAndSyst_ElePTsplit_2022preEE"]

###########################################
# Apply Scale Correction (for Data)       #
###########################################

# The scale correction is applied to data. Note: The corrections are deterministic.
# The evaluator returns a multiplicative factor that is applied to the pt using the example variables (run, ScEta, r9, abs(ScEta), pt, seedGain).
scale = scale_evaluator.evaluate("scale", data_run, data_electrons.ScEta, data_electrons.r9, np.abs(data_electrons.ScEta), data_electrons.pt, data_electrons.seedGain)
print(f"Multiplicative scale (for data): {scale}")

# Apply the scale to the electron pt
data_pt_corrected = scale * data_electrons.pt

#################################
# Apply Smearing Correction (MC)#
#################################

# For MC, we apply the smearing correction. The correction returns a parameter "smear" that indicates the needed smearing width.
# The model implemented in the smearing is a simple Gaussian smearing so that the effective correction factor per electron is 1 + N(0,1) * smear.

# Evaluate the nominal smearing width using example variables (pt, r9, and abs(ScEta))
smear = smear_and_syst_evaluator.evaluate("smear", mc_electrons.pt, mc_electrons.r9, np.abs(mc_electrons.ScEta))
print(f"Smearing width (for MC): {smear}")

# Calculate the nominal smearing factor.
# Since the smearing is stochastic, a random number is needed for each event.
random_numbers = rng.normal(loc=0.0, scale=1.0, size=len(mc_electrons.pt))
mc_pt_corrected_nominal = mc_electrons.pt * (1 + smear * random_numbers)

######################
# Uncertainties (MC) #
######################

# Obtain the uncertainty on the smearing width using MC original variables (pt, r9, and abs(ScEta)).
unc_smear = smear_and_syst_evaluator.evaluate("esmear", mc_electrons.pt, mc_electrons.r9, np.abs(mc_electrons.ScEta))
print(f"Smearing uncertainties: {unc_smear}")
mc_pt_corrected_smearing_up   = mc_electrons.pt * (1 + (smear + unc_smear) * random_numbers)
mc_pt_corrected_smearing_down = mc_electrons.pt * (1 + (smear - unc_smear) * random_numbers)
# Note: An alternative approach is to use the "smear_up" and "smear_down" variations, which are also provided for convenience.
# The relations are simple: smear_up = smear + unc_smear and smear_down = smear - unc_smear.

print("-------------------------------------------------------------------------------------")
print("--- Mean and Standard Deviation Comparison for MC electron pt (validate smearing) ---")
print(f"{'':<30} {'Original':<20} {'Nominal':<20} {'Smearing Up':<20} {'Smearing Down':<20}")
print(f"{'Mean':<30} {np.mean(ak.to_numpy(mc_electrons.pt)):<20.4f} {np.mean(mc_pt_corrected_nominal):<20.4f} {np.mean(mc_pt_corrected_smearing_up):<20.4f} {np.mean(mc_pt_corrected_smearing_down):<20.4f}")
print(f"{'Standard Deviation':<30} {np.std(ak.to_numpy(mc_electrons.pt)):<20.4f} {np.std(mc_pt_corrected_nominal):<20.4f} {np.std(mc_pt_corrected_smearing_up):<20.4f} {np.std(mc_pt_corrected_smearing_down):<20.4f}")
# Note that it behaves as expected: The mean is unchanged when smearing if the original distribution is Gaussian.
# However, the standard deviation increases when applying the smearing.
# The up/down variations sandwich the corrected (smeared) standard deviation.

# We now turn to the scale uncertainty, which is also evaluated on MC original variables (pt, r9, and abs(ScEta)) BUT applied on the smeared pt.
unc_scale = smear_and_syst_evaluator.evaluate("escale", mc_electrons.pt, mc_electrons.r9, np.abs(mc_electrons.ScEta))
mc_pt_corrected_scale_up   = (1 + unc_scale) * mc_pt_corrected_nominal
mc_pt_corrected_scale_down = (1 - unc_scale) * mc_pt_corrected_nominal
# Note: An alternative approach is to use the "scale_up" and "scale_down" variations, which are also provided for convenience.
# The relations are simple: scale_up = 1 + unc_scale and scale_down = 1 - unc_scale.

print("------------------------------------------------------------------------------------------------")
print("--- Mean and Standard Deviation Comparison for MC electron pt (validate scale uncertainties) ---")
print(f"{'':<30} {'Nominal (Smeared)':<20} {'Scale Up':<20} {'Scale Down':<20}")
print(f"{'Mean':<30} {np.mean(ak.to_numpy(mc_pt_corrected_nominal)):<20.4f} {np.mean(mc_pt_corrected_scale_up):<20.4f} {np.mean(mc_pt_corrected_scale_down):<20.4f}")
print(f"{'Standard Deviation':<30} {np.std(ak.to_numpy(mc_pt_corrected_nominal)):<20.4f} {np.std(mc_pt_corrected_scale_up):<20.4f} {np.std(mc_pt_corrected_scale_down):<20.4f}")
