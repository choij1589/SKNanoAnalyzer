#!/usr/bin/env python3
import os
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--era', type=str, default='2016preVFP', help='era')
args = parser.parse_args()

#This script is used to generate the path information for the sample data
#Sample information is stored in the CommonSampleInfo.json
#This script will generate the path information for the sample data and save it to the ForSNU folder
if args.era == "2016preVFP" or args.era == "2016postVFP":
    YEAR = "2016"
    BASEPATH = os.environ["SKNANO_RUN2_NANOAODPATH"]
elif args.era == "2017":
    YEAR = "2017"
    BASEPATH = os.environ["SKNANO_RUN2_NANOAODPATH"]
elif args.era == "2018":
    YEAR = "2018"
    BASEPATH = os.environ["SKNANO_RUN2_NANOAODPATH"]
elif args.era == "2022":
    YEAR = "2022"
    BASEPATH = os.environ["SKNANO_RUN3_NANOAODPATH"]
elif args.era == "2022EE":
    YEAR = "2022"
    BASEPATH = os.environ["SKNANO_RUN3_NANOAODPATH"]
elif args.era == "2023":
    YEAR = "2023"
    BASEPATH = os.environ["SKNANO_RUN3_NANOAODPATH"]
elif args.era == "2023BPix":
    YEAR = "2023"
    BASEPATH = os.environ["SKNANO_RUN3_NANOAODPATH"]
else:
    raise ValueError(f"Unknown era: {args.era}")

def parse_rootfiles_from(basePath):
    if not os.path.exists(basePath):
        print(f"No {basePath}")

    filePaths = []
    for root, _, files in os.walk(basePath):
        for file in files:
            if file.endswith(".root") and "wmLHE" not in file:
                filePaths.append(os.path.join(root, file))
    
    # sort the file paths by the number in the file name
    def extract_number(file_path):
        file_name = os.path.basename(file_path)
        # Remove the file extension
        file_name = os.path.splitext(file_name)[0]
        # Extract the number after the last underscore
        try:
            number = int(file_name.split('_')[-1])
            return number
        except ValueError:
            # If conversion fails, return the original string to maintain stable sorting
            return file_name

    filePaths.sort(key=extract_number)
    return filePaths

def main():
    sampleInfoJson = os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'CommonSampleInfo.json')
    sampleInfos = json.load(open(sampleInfoJson))
    os.makedirs(os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'ForSNU'), exist_ok=True)
    for alias, sampleInfo in sampleInfos.items():
        print(f"Processing {alias}...")
        # Check if the sample is MC or data
        isMC = sampleInfo["isMC"]
        if isMC:
            basePath = os.path.join(BASEPATH, args.era, sampleInfo["PD"])
            fileJsonPath = os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'ForSNU', alias+'.json')
            jsonData = {
                "name": alias,
                "isMC": 1,
                "PD": sampleInfo["PD"],
                "xsec": sampleInfo["xsec"],
                "sumsign": -1,
                "sumW": -1,
                "nmc": -1,
                "path": parse_rootfiles_from(basePath)
            }
            with open(fileJsonPath, 'w') as f:
                json.dump(jsonData, f, indent=4)
        else:
            for period in sampleInfo["periods"]:
                basePath = os.path.join(BASEPATH, args.era, f"{alias}/Run{YEAR}{period}")
                fileJsonPath = os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'ForSNU', alias+'_'+period+'.json')
                jsonData = {
                    "name": alias,
                    "isMC": 0,
                    "NEvents": [],  # Empty list for now, can be filled later
                    "path": parse_rootfiles_from(basePath)
                }
                with open(fileJsonPath, 'w') as f:
                    json.dump(jsonData, f, indent=4)


if __name__ == "__main__":
    main()
