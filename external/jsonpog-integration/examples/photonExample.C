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

    cout << key << " inputs:";
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input: sf->inputs()) {
        cout << ' ' << input.name() << flush;
        inputs.push_back(example.at(input.name()));
    }
    cout << endl;
    double result = sf->evaluate(inputs);
    cout << "JSON result: " << result << endl;
}

int main ()
{
    /*** example how to read the photon format v2 ***/

    fs::path fname = "../POG/EGM/2016postVFP_UL/photon.json.gz";
    cout << "Loading JSON file: " << fname << endl;
    assert(fs::exists(fname));
    unique_ptr<correction::CorrectionSet> cset =
                correction::CorrectionSet::from_file(fname.string());

    // general example
    map<string, correction::Variable::Type> example {
        {"year", "2016postVFP"},
        {"ValType", "sfup"},
        {"WorkingPoint", "Medium"},
        {"eta", 1.1},
        {"pt", 34.0},
        {"CSEVBin", "EBInc"},
        {"HasPixBin", "EBInc"},
    };
    run(cset, "UL-Photon-ID-SF"     , example);

    // CSEV vs PSV
    example["ValType"] = "sf";
    example["WorkingPoint"] = "Loose";
    run(cset, "UL-Photon-CSEV-SF"   , example);
    run(cset, "UL-Photon-PixVeto-SF", example);

    // upper and lower variations
    example["ValType"] = "sfup";
    run(cset, "UL-Photon-PixVeto-SF", example);
    example["ValType"] = "sfdown";
    run(cset, "UL-Photon-PixVeto-SF", example);

    /*** example how to read the photon format v3 ***/

    fname = "../POG/EGM/2023_Summer23/photon.json.gz";
    cout << "Loading JSON file: " << fname << endl;
    assert(fs::exists(fname));
    cset = correction::CorrectionSet::from_file(fname.string());
    example["year"] = "2023PromptC";

    // general example (nearly the same as above for 2016)
    example["phi"] = -1.8;
    example["ValType"] = "sfup";
    example["WorkingPoint"] = "Medium";
    run(cset, "Photon-ID-SF"     , example);

    // CSEV vs PSV
    example["ValType"] = "sf";
    example["eta"] = 1.2;
    example["R9"] = 0.85;
    example["WorkingPoint"] = "Loose";
    run(cset, "Photon-CSEV-SF"   , example);
    example["R9"] = 0.98;
    run(cset, "Photon-PixVeto-SF", example);

    // upper and lower variations
    example["ValType"] = "sfup";
    run(cset, "Photon-PixVeto-SF", example);
    example["ValType"] = "sfdown";
    run(cset, "Photon-PixVeto-SF", example);

    return 0;
}
