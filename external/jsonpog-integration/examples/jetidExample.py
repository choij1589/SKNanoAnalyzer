#! /usr/bin/env python
# Example of how to read the JetID JSON file (compatible with correctionlib 2.6.0+)
# For more information, see the README in
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
from correctionlib import _core

# Load CorrectionSet
fname = "../POG/JME/2022_Summer22/jetid.json.gz"
if fname.endswith(".json.gz"):
    import gzip
    with gzip.open(fname,'rt') as file:
        data = file.read().strip()
        evaluator = _core.CorrectionSet.from_string(data)
else:
    evaluator = _core.CorrectionSet.from_file(fname)


##### Jet Identification (JetID)
print(f"{' Barrel Jet ':=^80}")
eta = -0.23
chHEF = 0.4
neHEF = 0.25
chEmEF = 0.1
neEmEF = 0.15
muEF = 0.1
chMultiplicity = 1
neMultiplicity = 2
multiplicity = chMultiplicity + neMultiplicity
print(f"for a jet with eta: {eta},")
print(f"               chHEF: {chHEF}, neHEF: {neHEF}, chEmEF: {chEmEF}, neEmEF: {neEmEF}, muEF: {muEF},")
print(f"               chMultiplicity: {chMultiplicity}, neMultiplicity: {neMultiplicity}, multiplicity: {multiplicity}")

criteria_name = "AK4PUPPI_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4PUPPI_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

print("="*80)
print()

print(f"{' Barrel Jet with large muon energy fraction ':=^80}")
eta = -0.23
chHEF = 0.04
neHEF = 0.02
chEmEF = 0.02
neEmEF = 0.02
muEF = 0.9
chMultiplicity = 1
neMultiplicity = 2
multiplicity = chMultiplicity + neMultiplicity
print(f"for a jet with eta: {eta},")
print(f"               chHEF: {chHEF}, neHEF: {neHEF}, chEmEF: {chEmEF}, neEmEF: {neEmEF}, muEF: {muEF},")
print(f"               chMultiplicity: {chMultiplicity}, neMultiplicity: {neMultiplicity}, multiplicity: {multiplicity}")

criteria_name = "AK4PUPPI_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4PUPPI_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

print("="*80)
print()

print(f"{' Forward Jet ':=^80}")
eta = 4.5
chHEF = 0.4
neHEF = 0.25
chEmEF = 0.1
neEmEF = 0.15
muEF = 0.1
chMultiplicity = 1
neMultiplicity = 2
multiplicity = chMultiplicity + neMultiplicity
print(f"for a jet with eta: {eta},")
print(f"               chHEF: {chHEF}, neHEF: {neHEF}, chEmEF: {chEmEF}, neEmEF: {neEmEF}, muEF: {muEF},")
print(f"               chMultiplicity: {chMultiplicity}, neMultiplicity: {neMultiplicity}, multiplicity: {multiplicity}")

criteria_name = "AK4PUPPI_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4PUPPI_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_Tight"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

print("="*80)
print()

