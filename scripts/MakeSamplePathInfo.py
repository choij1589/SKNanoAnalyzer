import os
import json

#This script is used to generate the path information for the sample data
#Sample information is stored in the CommonSampleInfo.json
#This script will generate the path information for the sample data
#basePath = '/gv0/Users/yeonjoon/DATA/SKFlat/Run3NanoAODv12/' #currently my private directory
basePath = "/gv0/DATA/SKNano/Run3NanoAODv12"
eras = ['2022','2022EE']
for era in eras:
    sampleInfoJson = os.path.join(os.environ['SKNANO_DATA'],era,'Sample','CommonSampleInfo.json')
    sampleInfos = json.load(open(sampleInfoJson))
    isMC = 'MC'
    for alias,sampleInfo in sampleInfos.items():
        path = os.path.join(basePath,era,isMC,sampleInfo['PD'])
        filePaths = []
        #Folder structure is not fixed yet, so let's do the recursive search until the .root file appears, and save all absolute paths
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith('.root'):
                    filePaths.append(os.path.join(root,file))
        
        #now save the path information to another json file
        fileJsonPath = os.path.join(os.environ['SKNANO_DATA'],era,'Sample','ForSNU',alias+'.json')
        with open(fileJsonPath,'w') as f:
            json.dump(filePaths, f, indent=4)
