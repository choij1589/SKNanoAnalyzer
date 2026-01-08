#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include "correction.h"

using namespace std;
namespace fs = filesystem;

void run (const unique_ptr<correction::CorrectionSet>& cset,
          const string& key,
          const map<string, correction::Variable::Type>& example)
{
    correction::Correction::Ref sf = cset->at(key);
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input: sf->inputs())
        inputs.push_back(example.at(input.name()));
    double result = sf->evaluate(inputs);
    cout << "JSON result: " << result << endl;
}

int main ()
{
    ////////////////////////////
    // Example A: 2016postVFP //
    ////////////////////////////

    fs::path fname = "../POG/MUO/2016postVFP_UL/muon_Z.json.gz";
    cout << "Loading JSON file: " << fname << endl;
    assert(fs::exists(fname));
    unique_ptr<correction::CorrectionSet> cset =
                correction::CorrectionSet::from_file(fname.string());

    // TrackerMuon Reconstruction UL scale factor ==> NOTE the year key has been removed, for consistency with Run 3
    map<string, correction::Variable::Type> example {
        {"pt", 50.0}, // muon transverse momentum
        {"eta", -1.1}, // muon absolute pseudorapidity
        {"scale_factors", "nominal"}, // variation
    };
    run(cset, "NUM_TrackerMuons_DEN_genTracks", example);

    // Medium ID UL scale factor, down/up variations ==> NOTE the year key has been removed, for consistency with Run 3
    example = {{"eta", 0.8}, {"pt", 35.0}, {"scale_factors", "systdown"}};
    run(cset, "NUM_MediumID_DEN_TrackerMuons", example);
    example["scale_factors"] = "systup";
    run(cset, "NUM_MediumID_DEN_TrackerMuons", example);

    // Trigger UL systematic uncertainty only ==> NOTE the year key has been removed, for consistency with Run 3
    example = {{"eta", 1.8}, {"pt", 54.0}, {"scale_factors", "syst"}};
    run(cset, "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight", example);

    //////////////////////////
    // Example B: 2022preEE //
    //////////////////////////

    fname = "../POG/MUO/2022_Summer22/muon_Z.json.gz";
    cout << "Loading JSON file: " << fname << endl;
    assert(fs::exists(fname));
    cset = correction::CorrectionSet::from_file(fname.string());

    // Medium ID 2022 scale factor using eta as input
    example = {{"eta", -1.1}, {"pt", 45.0}, {"scale_factors", "nominal"}};
    run(cset, "NUM_MediumID_DEN_TrackerMuons", example);

    // Medium ID 2022 scale factor using eta as input ==> Note that this value should be the same
    // as the previous one, since even though the input can be signed eta, the SFs for 2022 were
    // computed for |eta|. This is valid for ALL the years and jsons
    example["eta"] = 1.1;
    run(cset, "NUM_MediumID_DEN_TrackerMuons", example);

    // Trigger 2022 systematic uncertainty only 
    example = {{"eta", -1.8}, {"pt", 54.0}, {"scale_factors", "syst"}};
    run(cset, "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", example);

    return 0;
}
