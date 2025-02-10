import ROOT
import correctionlib.convert as convert
import json
import argparse
import os
import concurrent.futures
from tqdm import tqdm

# ---------------------------------------------------------------------
# 1. Parse histogram keys from a ROOT file
# ---------------------------------------------------------------------
def histParser(root_file):
    """
    Parse all histogram keys in the given ROOT file.
    Expected key format (example):
      tagging#b##era#2022EE##tagger#deepJet##working_point#M##flavor#5##systematic#central##num
      or for the denominator:
      tagging#b##era#2022EE##flavor#5##systematic#central##den
    Returns a list of dictionaries with extracted fields and the full key.
    """
    hist_list = []
    for keyObj in root_file.GetListOfKeys():
        key = keyObj.GetName()
        if "##" not in key:
            continue
        segments = key.split("##")
        info = {}
        # Process each segment: if it contains '#' then split into field and value.
        for seg in segments:
            if '#' not in seg:
                # Fallback: sometimes the type might be appended at the end.
                info["hist_type"] = seg.strip()
                continue
            parts = seg.split("#", 1)
            if len(parts) < 2:
                continue
            field, value = parts
            info[field] = value
        # Try to set hist_type if not set: look at the key ending.
        if "hist_type" not in info:
            if key.endswith("num"):
                info["hist_type"] = "num"
            elif key.endswith("den"):
                info["hist_type"] = "den"
            elif key.endswith("eff"):
                info["hist_type"] = "eff"
        info["hist_key"] = key
        hist_list.append(info)
    return hist_list

# ---------------------------------------------------------------------
# 2. Create efficiency histograms by matching numerator and denominator histograms
# ---------------------------------------------------------------------
def makeTempEffHist(in_file):
    """
    From the input ROOT file (in_file), find matching numerator and denominator histograms.
    
    Note that the numerator histogram key contains extra fields (tagger and working_point)
    while the denominator is shared among taggers/working points. We match based on:
    
      - tagging
      - era
      - flavor
      - systematic
    
    For each numerator histogram, the corresponding denominator is the one whose key
    matches these fields and has hist_type "den".
    
    The resulting efficiency histogram (with key "eff") is saved to "output.root".
    """
    hist_list = histParser(in_file)
    
    # Separate numerator and denominator histograms
    num_list = [h for h in hist_list if h.get("hist_type", "") == "num"]
    den_list = [h for h in hist_list if h.get("hist_type", "") == "den"]

    eff_hists = []
    # Matching criteria for numerator vs. denominator:
    # Only match on these fields: tagging, era, flavor, systematic.
    match_fields = ["tagging", "era", "flavor", "systematic"]
    for num_info in num_list:
        # For each numerator, try to find a denominator that matches on the common fields.
        matching_den = None
        for den_info in den_list:
            print(f"Will find denominator for numerator: {[num_info[f] for f in match_fields]}")
            print("Total list of denominator: ")
            for den in den_list:
                print(f"Denominator: {[den[f] for f in match_fields]}")
            if all(num_info.get(f) == den_info.get(f) for f in match_fields):
                matching_den = den_info
                break
        if not matching_den:
            # If no matching denominator is found, skip this numerator.
            print(f"WARNING: No denominator found for numerator: {num_info['hist_key']}")
            continue

        num_hist = in_file.Get(num_info["hist_key"])
        den_hist = in_file.Get(matching_den["hist_key"])
        if not num_hist or not den_hist:
            print(f"WARNING: Missing histogram object for {num_info['hist_key']} or {matching_den['hist_key']}")
            continue

        # Create efficiency histogram by cloning numerator and dividing by denominator.
        eff_name = num_info["hist_key"].replace("num", "eff")
        eff_hist = num_hist.Clone(eff_name)
        eff_hist.Divide(den_hist)
        eff_hists.append(eff_hist)
        print(f"Created efficiency histogram: {eff_name}")

    # Write all efficiency histograms to a new ROOT file "output.root"
    out_file = ROOT.TFile("output.root", "RECREATE")
    out_file.cd()
    for hist in eff_hists:
        hist.Write()
    out_file.Close()

# ---------------------------------------------------------------------
# 3. Worker function to convert an efficiency histogram to correction data
# ---------------------------------------------------------------------
def process_histogram(hist_info):
    """
    Reopen the "output.root" file and convert the efficiency histogram given by hist_info.
    Returns a tuple: (tagger, systematic, working_point, flavor, conversion_data)
    """
    hist_path = f"output.root:{hist_info['hist_key']}"
    try:
        conv = convert.from_uproot_THx(hist_path)
        conv_data = conv.dict()['data']
        # Adjust conversion dictionary as in your original Code 2.
        conv_data['inputs'] = ['abseta', 'pt']
        conv_data['flow'] = 'clamp'
    except Exception as e:
        raise RuntimeError(f"Conversion failed for {hist_info['hist_key']}: {e}")

    # For grouping, extract tagger, working_point, and flavor.
    # (For the denominator these fields are missing, but we only process efficiency histograms,
    #  which come from numerators and thus include tagger and working_point.)
    tagger = hist_info.get("tagger")
    systematic = hist_info.get("systematic")
    working_point = hist_info.get("working_point")
    try:
        flavor = int(hist_info.get("flavor"))
    except (TypeError, ValueError):
        flavor = None

    return (tagger, systematic, working_point, flavor, conv_data)

# ---------------------------------------------------------------------
# 4. Build the JSON corrections file dynamically
# ---------------------------------------------------------------------
def makingJson(era, tagging_mode):
    """
    Open the merged efficiency ROOT file ("output.root"), parse its histogram keys,
    filter for a given era and tagging mode ("b" or "c"), convert each histogram in parallel,
    and group the results into a nested JSON structure.
    """
    out_file = ROOT.TFile("output.root", "READ")
    hist_list = histParser(out_file)
    out_file.Close()

    # Filter for efficiency histograms matching the era and tagging mode.
    filtered = [h for h in hist_list 
                if h.get("hist_type", "") == "eff" and
                   h.get("era") == era and
                   h.get("tagging") == tagging_mode]

    conversion_results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_histogram, h) for h in filtered]
        for future in tqdm(concurrent.futures.as_completed(futures),
                           total=len(futures),
                           desc="Converting histograms"):
            try:
                conversion_results.append(future.result())
            except Exception as e:
                print(f"Error converting histogram: {e}")

    # Group the conversion results by tagger → systematic → working_point → flavor.
    grouped = {}
    for tagger, systematic, wp, flavor, data in conversion_results:
        if None in (tagger, systematic, wp, flavor):
            continue
        grouped.setdefault(tagger, {}).setdefault(systematic, {}).setdefault(wp, {})[flavor] = data

    corrections = []
    for tagger, syst_dict in grouped.items():
        corr_entry = {
            "name": tagger,
            "version": 0,
            "inputs": [
                {"name": "systematic", "type": "string"},
                {"name": "working_point", "type": "string", "description": "b-tagging working point"},
                {"name": "flavor", "type": "int", "description": "hadron flavor definition: 5=b, 4=c, 0=udsg"},
                {"name": "abseta", "type": "real"},
                {"name": "pt", "type": "real"}
            ],
            "output": {"name": "eff", "type": "real"},
            "data": {
                "nodetype": "category",
                "input": "systematic",
                "content": []
            }
        }
        for systematic, wp_dict in syst_dict.items():
            syst_entry = {
                "key": systematic,
                "value": {
                    "nodetype": "category",
                    "input": "working_point",
                    "content": []
                }
            }
            for wp, flav_dict in wp_dict.items():
                wp_entry = {
                    "key": wp,
                    "value": {
                        "nodetype": "category",
                        "input": "flavor",
                        "content": []
                    }
                }
                for flav, data in flav_dict.items():
                    wp_entry["value"]["content"].append({"key": flav, "value": data})
                wp_entry["value"]["content"].sort(key=lambda x: x["key"])
                syst_entry["value"]["content"].append(wp_entry)
            syst_entry["value"]["content"].sort(key=lambda x: x["key"])
            corr_entry["data"]["content"].append(syst_entry)
        corr_entry["data"]["content"].sort(key=lambda x: x["key"])
        corrections.append(corr_entry)

    main_json = {
        "schema_version": 2,
        "description": "This json file contains the b-tagging efficiency corrections",
        "corrections": corrections
    }
    return main_json

# ---------------------------------------------------------------------
# 5. Main script: process the input file and write JSON files
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Create b-tagging efficiency JSON file dynamically")
    parser.add_argument('--input_folder', dest='input_folder',
                        help="Folder containing ROOT files",
                        default=os.path.join(os.environ.get('SKNANO_OUTPUT', ''), 'MeasureJetTaggingEff'))
    parser.add_argument('--input', dest='input',
                        help="Input ROOT file name", default='TTLJ_powheg.root')
    parser.add_argument('--out_name_str', dest='out_name_str',
                        help="Output JSON file name prefix", default='')
    args = parser.parse_args()

    # List of eras to process (adjust as needed)
    totalEras = ['2017']
    for era in totalEras:
        out_dir = os.path.join(os.environ.get('SKNANO_DATA', ''), era, 'BTV')
        os.makedirs(out_dir, exist_ok=True)
        input_path = os.path.join(args.input_folder, era, args.input)
        in_file = ROOT.TFile(input_path)
        if not in_file or in_file.IsZombie():
            raise RuntimeError(f"Cannot open file: {input_path}")

        # Create efficiency histograms by matching numerator and denominator histograms.
        makeTempEffHist(in_file)
        in_file.Close()

        # Process both tagging modes: "b" and "c"
        for tagging_mode, json_suffix in zip(["b", "c"], ["btaggingEff.json", "ctaggingEff.json"]):
            main_json = makingJson(era, tagging_mode)
            out_path = os.path.join(out_dir, args.out_name_str + json_suffix)
            with open(out_path, "w") as fout:
                json.dump(main_json, fout, indent=4)
            print(f"Wrote JSON file: {out_path}")

if __name__ == "__main__":
    main()
