#!/usr/bin/env python3
import json
import os
import argparse
import ROOT
from pprint import pprint

parser = argparse.ArgumentParser()
parser.add_argument("--era", type=str, required=True, help="era")
args = parser.parse_args()
SKNanoOutputDir = os.environ["SKNANO_OUTPUT"]
SKNanoDataDir = os.environ["SKNANO_DATA"]

with open(f"{SKNanoDataDir}/{args.era}/Sample/CommonSampleInfo.json", "r") as f:
    common_info = json.load(f)

# Here we assume that we already parsed the sample paths using scripts/MakeSamplePathInfo.py
# and run GetEffLumi for all samples to get the inputs for this script.

# For data samples, only have to parse the number of events
def parseDataInfoFor(sample_name):
    # parse periods from sample_info
    periods = common_info[sample_name]["periods"]
    nevts = []
    for period in periods:
        try:
            f = ROOT.TFile.Open(f"{SKNanoOutputDir}/GetEffLumi/{args.era}/{sample_name}_{period}.root")
            h = f.Get("NEvents")
            nevts.append(h.GetBinContent(1))
            f.Close()
        except:
            print(f"Error opening file for {sample_name} {period}")
            nevts.append(-1)
    
    common_info[sample_name]["NEvents"] = nevts

    # add in ForSNU/$SAMPLE_NAME_PERIOD.json
    for period in periods:
        with open(f"{SKNanoDataDir}/{args.era}/Sample/ForSNU/{sample_name}_{period}.json", "r") as f:
            for_snu_info = json.load(f)
        for_snu_info["NEvents"] = nevts
        with open(f"{SKNanoDataDir}/{args.era}/Sample/ForSNU/{sample_name}_{period}.json", "w") as f:
            json.dump(for_snu_info, f, indent=4)

def parseMCInfoFor(sample_name):
    nevts = 0
    sumsign = 0
    sumW = 0
    filterEff = None
    try:
        f = ROOT.TFile.Open(f"{SKNanoOutputDir}/GetEffLumi/{args.era}/{sample_name}.root")
        h = f.Get("NEvents")
        nevts = h.GetBinContent(1)
        h = f.Get("sumSign")
        sumsign = h.GetBinContent(1)
        h = f.Get("sumW")
        sumW = h.GetBinContent(1)
        h_passed = f.Get("GenFilter_numEventsPassed")
        h_total = f.Get("GenFilter_numEventsTotal")
        if h_passed and h_total:
            totalPassed = h_passed.GetBinContent(1)
            totalTotal = h_total.GetBinContent(1)
            if totalTotal > 0:
                filterEff = totalPassed / totalTotal
        f.Close()
    except:
        print(f"Error opening file for {sample_name}")
        sumsign = -1
        sumW = -1
    common_info[sample_name]["nmc"] = nevts
    common_info[sample_name]["sumsign"] = sumsign
    common_info[sample_name]["sumW"] = sumW
    if filterEff is not None:
        common_info[sample_name]["filterEff"] = filterEff

    # add in ForSNU/$SAMPLE_NAME.json
    with open(f"{SKNanoDataDir}/{args.era}/Sample/ForSNU/{sample_name}.json", "r") as f:
        for_snu_info = json.load(f)
    for_snu_info["nmc"] = nevts
    for_snu_info["sumsign"] = sumsign
    for_snu_info["sumW"] = sumW
    with open(f"{SKNanoDataDir}/{args.era}/Sample/ForSNU/{sample_name}.json", "w") as f:
        json.dump(for_snu_info, f, indent=4)

def main():
    for sample_name in common_info.keys():
        if common_info[sample_name]["isMC"]:
            parseMCInfoFor(sample_name)
        else:
            parseDataInfoFor(sample_name)
    
    # dump the updated sample_info
    with open(f"{SKNanoDataDir}/{args.era}/Sample/CommonSampleInfo.json", "w") as f:
        json.dump(common_info, f, indent=4)

if __name__ == "__main__":
    main()
