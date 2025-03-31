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
BASEPATH = '/gv0/Users/choij/SKNano/' 
if args.era == "2016preVFP" or args.era == "2016postVFP":
    YEAR = "2016"
elif args.era == "2017":
    YEAR = "2017"
elif args.era == "2018":
    YEAR = "2018"

def parse_rootfiles_from(basePath):
    filePaths = []
    for root, _, files in os.walk(basePath):
        for file in files:
            if file.endswith(".root"):
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
    for alias, sampleInfo in sampleInfos.items():
        print(f"Processing {alias}...")
        # Check if the sample is MC or data
        isMC = sampleInfo["isMC"]
        if isMC:
            basePath = os.path.join(BASEPATH, sampleInfo["PD"], f"crab_MC_{args.era}_{sampleInfo['PD']}")
            fileJsonPath = os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'ForSNU', alias+'.json')
            jsonData = {
                "name": alias,
                "isMC": 1,
                "PD": sampleInfo["PD"],
                "xsec": sampleInfo["xsec"],
                "sumsign": -1,
                "sumW": -1,
                "path": parse_rootfiles_from(basePath)
            }
            with open(fileJsonPath, 'w') as f:
                json.dump(jsonData, f, indent=4)
        else:
            for period in sampleInfo["periods"]:
                basePath = os.path.join(BASEPATH, alias, f"crab_DATA_{args.era}_{alias}_Run{YEAR}{period}")
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