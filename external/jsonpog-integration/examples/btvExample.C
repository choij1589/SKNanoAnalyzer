#include <cassert>
#include <cstdlib>
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

double rndm (double scale = 1, double m = 0)
{
    return m + (scale * rand() ) / RAND_MAX;
}

int GetFlId (double n)
{
    if      (n > 0.8) return 5; // b jet
    else if (n > 0.6) return 4; // c jet
    else              return 0; // light jet
}

int main ()
{
    fs::path fname = "../POG/BTV/2018_UL/btagging.json.gz";
    cout << "Loading JSON file: " << fname << endl;
    assert(fs::exists(fname));
    unique_ptr<correction::CorrectionSet> cset =
                correction::CorrectionSet::from_file(fname.string());

    // generate 20 dummy jet features
    for (int i = 0; i < 20; ++i) {

        double pt = rndm(900, 20),
               abseta = rndm(2.4),
               discriminant = rndm();
        int flavor = GetFlId(rndm());
        cout << "========\n" << pt << ' ' << abseta << ' ' << discriminant << ' ' << flavor << endl;

        map<string, correction::Variable::Type> example {
            /* jet properties */
            {"pt"  , pt}, // jet transverse momentum
            {"abseta" , abseta}, // absolute jet pseudorapidity
            {"flavor", flavor}, // jet flavour
            {"discriminant", discriminant}, // jet discriminant
            /* analysis dependent */
            {"systematic", "central"}, // systematic variation
            {"working_point", "M"}, // discriminant working point
        };

        /*** case 1: fixedWP correction with mujets (here medium WP) ***/
        if (flavor != 0) run(cset, "deepJet_mujets", example);
        else               run(cset, "deepJet_incl"  , example);

        /*** case 2: fixedWP correction uncertainty (here tight WP and comb SF) ***/
        example["systematic"] = "up_correlated";
        example["working_point"] = "T";
        if (flavor != 0) run(cset, "deepJet_comb", example);
        else               run(cset, "deepJet_incl", example);

        /*** case 3: shape correction SF ***/
        example["systematic"] = "central";
        run(cset, "deepJet_shape", example);

        /*** case 4: shape correction SF uncertainties ***/
        if (flavor != 4) {
            example["systematic"] = "up_hfstats2";
            run(cset, "deepJet_shape", example);
        }
        else {
            example["systematic"] = "up_cferr1";
            run(cset, "deepJet_shape", example);
        }
    }

    return 0;
}
