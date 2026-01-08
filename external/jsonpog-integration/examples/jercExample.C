#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include "correction.h"

using namespace std;
namespace fs = filesystem;

double singleLevel (const map<string, correction::Variable::Type>& example,
                    const unique_ptr<correction::CorrectionSet>& cset,
                    string jec, string lvl, string algo)
{
    string key = jec + '_' + lvl + '_' + algo;
    cout << "JSON access to key: " << key << endl;

    correction::Correction::Ref sf = cset->at(key);

    cout << "Inputs:";
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input: sf->inputs()) {
        cout << ' ' << input.name();
        inputs.push_back(example.at(input.name()));
    }
    cout << endl;

    double result = sf->evaluate(inputs);
    cout << "JSON result: " << result << endl;

    return result;
}

double compoundLevel (const map<string, correction::Variable::Type>& example,
                      const unique_ptr<correction::CorrectionSet>& cset,
                      string jec, string lvl, string algo)
{
    string key = jec + '_' + lvl + '_' + algo;
    cout << "JSON access to key: " << key << endl;

    correction::CompoundCorrection::Ref sf = cset->compound().at(key); // note: the only different is here

    cout << "Inputs:";
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input: sf->inputs()) {
        cout << ' ' << input.name();
        inputs.push_back(example.at(input.name()));
    }
    cout << endl;

    double result = sf->evaluate(inputs);
    cout << "JSON result: " << result << endl;

    return result;
}

void smear (const map<string, correction::Variable::Type>& example,
            const unique_ptr<correction::CorrectionSet>& cset)
{
    correction::Correction::Ref sf = cset->at("JERSmear");

    cout << "Inputs:";
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input: sf->inputs()) {
        cout << ' ' << input.name();
        inputs.push_back(example.at(input.name()));
    }
    cout << endl;

    double result = sf->evaluate(inputs);
    cout << "JSON result: " << result << endl;

    // to implement smearing in the analysis code, multiply
    // the `jersmear_factor` obtained above to the `Jet_pt`
    // and `Jet_mass` variables
}

int main ()
{
    map<string, correction::Variable::Type> example {
        { // jet transverse momentum
        "JetPt", 100.0},
        { // jet pseudorapidity
        "JetEta", 0.0},
        { // jet azimuthal angle
        "JetPhi", 0.2},
        { // jet area
        "JetA", 0.5},
        { // median energy density (pileup)
        "Rho", 15.0},
        { // systematic variation (only for JER SF)
        "systematic", "nom"},
        { // pT of matched gen-level jet (only for JER smearing)
        "GenPt", 80.0}, // or -1 if no match
        { // unique event ID used for deterministic pseudorandom number generation (only for JER smearing)
        "EventID", 12345},
    };

    string jec = "Summer19UL16_V7_MC", // JEC base tag
           jer = "Summer20UL16_JRV3_MC", // JER base tag
           algo_ak4 = "AK4PFchs", // AK4 jet algorithm
           algo_ak8 = "AK8PFPuppi", // AK8 jet algorithm
           lvl_single = "L2Relative", // jet energy correction level
           lvl_compound = "L1L2L3Res", // jet energy correction level
           unc = "Total"; // jet energy uncertainty

    // print input information
    cout << "\nJEC parameters\n##############"
         << "\njec = " << jec
         << "\nalgo_ak4 = " << algo_ak4
         << "\nalgo_ak8 = " << algo_ak8;
    //cout << "\nJetPt = " << get<double>(example["JetPt"]) << endl;
    for (const string& v: {"JetPt", "JetEta", "JetA", "JetPhi", "JetA", "Rho"})
        cout << '\n' << v << ' ' << get<double>(example[v]);
    cout << '\n' << endl;

    /**** load JSON files using correctionlib ****/

    // AK4
    fs::path fname_ak4 = "../POG/JME/2016postVFP_UL/jet_jerc.json.gz";
    cout << "Loading JSON file: " << fname_ak4 << endl;
    assert(fs::exists(fname_ak4));
    unique_ptr<correction::CorrectionSet> cset_ak4 =
                correction::CorrectionSet::from_file(fname_ak4.string());

    // AK8
    fs::path fname_ak8 = "../POG/JME/2016postVFP_UL/fatJet_jerc.json.gz";
    cout << "Loading JSON file: " << fname_ak8 << endl;
    assert(fs::exists(fname_ak8));
    unique_ptr<correction::CorrectionSet> cset_ak8 =
                correction::CorrectionSet::from_file(fname_ak8.string());

    // tool for JER smearing
    fs::path fname_jersmear = "../POG/JME/jer_smear.json.gz";
    cout << "Loading JSON file: " << fname_jersmear << endl;
    assert(fs::exists(fname_jersmear));
    auto cset_jersmear = correction::CorrectionSet::from_file(fname_jersmear.string());

    /**** run examples ****/

    cout << "\nExample 1: single JEC level\n===================" << endl;
    singleLevel(example, cset_ak4, jec, lvl_single, algo_ak4);
    singleLevel(example, cset_ak8, jec, lvl_single, algo_ak8);

    cout << "\nExample 2: compound JEC level\n===================" << endl;
    compoundLevel(example, cset_ak4, jec, lvl_compound, algo_ak4);
    compoundLevel(example, cset_ak8, jec, lvl_compound, algo_ak8);

    cout << "\nExample 3: JEC uncertainty source\n===================" << endl;
    singleLevel(example, cset_ak4, jec, unc, algo_ak4);
    singleLevel(example, cset_ak8, jec, unc, algo_ak8);
    // additional note: Regrouped/reduced set of uncertainty sorces as detailed in
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Run_2_reduced_set_of_uncertainty  # noqa
    // are included in relevant JSON files (currently UL) with a "Regrouped_"-prefix,
    // e.g. for 2016 one could access "Absolute_2016" via:
    // sf = cset["Summer19UL16_V7_MC_Regrouped_Absolute_2016_AK4PFchs"]

    cout << "\nExample 4: JER scale factor\n===================" << endl;
    example["JERSF"] = singleLevel(example, cset_ak4, jer, "ScaleFactor", algo_ak4);

    cout << "\nExample 5: JER (pT resolution)\n===================" << endl;
    example["JER"] = singleLevel(example, cset_ak4, jer, "PtResolution", algo_ak4);

    cout << "\nExample 6: JER smearing\n===================" << endl;
    smear(example, cset_jersmear);

    return 0;
}
