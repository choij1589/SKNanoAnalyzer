#!/usr/bin/env python3
#This it the preliminary version of SKFlat.py
#Using htcondor python binding and DAGMAN workflow manager
import os, shutil
import warnings
import argparse
import htcondor
import datetime
import json
import re

from htcondor import dags
from tqdm.rich import tqdm
from tqdm import TqdmExperimentalWarning

from templates.job_dict import main_job, hadd_job, final_job

dag = dags.DAG()
warnings.simplefilter("ignore", TqdmExperimentalWarning)

##############################
#Global Variables
##############################
SKNANO_HOME = os.environ['SKNANO_HOME']
SKNANO_RUNLOG = os.environ['SKNANO_RUNLOG']
SKNANO_OUTPUT = os.environ['SKNANO_OUTPUT']
SKNANO_DATA = os.environ['SKNANO_DATA']
SKNANO_LIB = os.environ['SKNANO_LIB']
SKNANO_INSTALLDIR = os.environ['SKNANO_INSTALLDIR']
SKNANO_RUN3_NANOAODPATH = os.environ['SKNANO_RUN3_NANOAODPATH']
SKNANO_RUN2_NANOAODPATH = os.environ['SKNANO_RUN2_NANOAODPATH']
username = os.environ['USER']
Run = {'2016preVFP':2,'2016postVFP':2,'2017':2,'2018':2,'2022':3,'2022EE':3, '2023':3, '2023BPix':3}
TOKEN = os.environ['TOKEN_TELEGRAMBOT']
chat_id = os.environ['USER_CHATID']
url = f"https://api.telegram.org/bot{TOKEN}/sendMessage?chat_id={chat_id}"
SKIMMING_MODE = False

##############################
#Load commonSampleInfo.json at start
sampleInfoJsons = {}
for era in Run.keys():
    sampleInfoJsons[era] = json.load(open(os.path.join(SKNANO_DATA,era,'Sample','CommonSampleInfo.json')))
skimInfoJsons = {}


for era in Run.keys():
    try:
        skimInfoJsons[era] = json.load(open(os.path.join(SKNANO_DATA,era,'Sample','Skim','skimTreeInfo.json')))
    except:
        print(f"\033[93mWarning: {era} skimTreeInfo.json is not exist\033[0m")

##############################
def isMCandGetPeriod(sample):
    #if sample is ends with _one capital letter, it is data
    #hope there will be no exception(please)
    if sample.split("_")[-1].isupper() and len(sample.split("_")[-1]) == 1:
        return False, sample.split("_")[-1]
    else:
        return True, None

def getSkimmingOutBaseAndSuffix(era, sample, AnalyzerName):
    isMC, period = isMCandGetPeriod(sample)
    suffix = f"Temp_Skim_{AnalyzerName.replace('Skim_','')}_{sample if isMC else sample.replace(f'_{period}','')}"
    if Run[era] == 2:
        out_base = os.path.join(SKNANO_RUN2_NANOAODPATH ,era,'MC' if isMC else 'DATA','Skim',os.environ['USER'],suffix,'' if isMC else f'Period{period}', 'tree.root') if SKIMMING_MODE else 'output/hists.root'
    elif Run[era] == 3:
        out_base = os.path.join(SKNANO_RUN3_NANOAODPATH ,era,'MC' if isMC else 'DATA','Skim',os.environ['USER'],suffix,'' if isMC else f'Period{period}', 'tree.root') if SKIMMING_MODE else 'output/hists.root'

    return out_base, suffix

def getEraList(eras, runs):
    if runs == 'None':
        if eras == 'All':
            return Run.keys() 
        eras = eras.split(",")
        for era in eras:
            if era not in Run.keys():
                print('\033[91m'+f"ERROR: {era} is not a valid era"+'\033[0m')
                exit()
    else:
        print('\033[93m'+"Warning: --era option will be ignored because --run option is given"+'\033[0m')
        runs = runs.split(",")
        runs = [run.replace(" ","") for run in runs]
        eras = []
        if 'Run2' in runs:
            eras += [e for e, r in Run.items() if r == 2]
        if 'Run3' in runs:
            eras += [e for e, r in Run.items() if r == 3]
    return eras

def makeSampleList(samplelist,era):
    #add Period to data sample, and wildcard search
    copylist = []
    for sample in samplelist:
        if sample.startswith("Skim_"):
            if SKIMMING_MODE:
                print('\033[91m'+f"Error: Skimmed sample {sample} is not allowed in skimming mode"+'\033[0m')
                exit()
            if '*' in sample:
                #wildcard search using regex
                sample = sample.replace('*','.*')
                sample = re.compile(sample)
                for sampleInfo in skimInfoJsons[era].keys():
                    if sample.match(sampleInfo):
                        if skimInfoJsons[era][sampleInfo]['isMC']:
                            copylist.append(sampleInfo)
                        else:
                            copylist += [f"{sampleInfo}_{period}" for period in skimInfoJsons[era][sampleInfo]['periods']]
                continue
            elif sample not in skimInfoJsons[era]:
                print('\033[93m'+f"Warning: {sample} is not exist in era {era}"+'\033[0m')
                continue
            skimInfo = skimInfoJsons[era][sample]
            sampleInfo = sampleInfoJsons[era][skimInfo['PD']]
            if sampleInfo['isMC']:
                copylist.append(sample)
            if not sampleInfo['isMC']:
                copylist += [f"{sample}_{period}" for period in sampleInfo['periods']]
        else:
            if '*' in sample:
                #wildcard search using regex
                sample = sample.replace('*','.*')
                sample = re.compile(sample)
                for sampleInfo in sampleInfoJsons[era].keys():
                    if sample.match(sampleInfo):
                        if sampleInfoJsons[era][sampleInfo]['isMC']:
                            copylist.append(sampleInfo)
                        else:
                            copylist += [f"{sampleInfo}_{period}" for period in sampleInfoJsons[era][sampleInfo]['periods']]
                continue
            elif sample not in sampleInfoJsons[era]:
                print('\033[93m'+f"Warning: {sample} is not exist in era {era}"+'\033[0m')
                continue
            sampleInfo = sampleInfoJsons[era][sample]
            if sampleInfo['isMC']:
                copylist.append(sample)
            if not sampleInfo['isMC']:
                copylist += [f"{sample}_{period}" for period in sampleInfo['periods']]

    return copylist

def getUserFlagsList(Userflags):
    UserflagsList = Userflags.split(",")
    UserflagsList = [x for x in UserflagsList if x != ""]
    return UserflagsList

def getTimeStamp():
    ## TimeStamp

    # 1) dir/file name style
    JobStartTime = datetime.datetime.now()
    timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')
    timestring =  JobStartTime.strftime('%Y-%m-%d %H:%M:%S')
    return timestamp, timestring

def setParser():
    parser = argparse.ArgumentParser(description='SKNano Command')
    parser.add_argument('-a', dest='Analyzer', default="")
    parser.add_argument('-i', dest='InputSample', default="", help="Input sample list, can be txt file, or comma separated")
    #parser.add_argument('-p', dest='DataPeriod', default="ALL")
    #parser.add_argument('-l', dest='InputSampleList', default="")
    parser.add_argument('-n', dest='NJobs', default=1, type=int, help="Number of jobs, If Negative, it will means number of files per job")
    #parser.add_argument('-o', dest='Outputdir', default="")
    #parser.add_argument('-q', dest='Queue', default="fastq")
    parser.add_argument('-e', dest='Era', default="All",help="2022, 2022EE. can be comma separated")
    parser.add_argument('-r', dest='Run', default="None",help="Run2, Run3. can be comma separated. override era option")
    parser.add_argument('--userflags', dest='Userflags', default="")
    parser.add_argument('--nmax', dest='NMax', default=300, type=int, help="maximum running jobs")
    parser.add_argument('--reduction', dest='Reduction', default=1, type=float)
    parser.add_argument('--memory', dest='Memory', default=2048, type=float)
    parser.add_argument('--ncpu', dest='ncpu', default=1, type=int) 
    parser.add_argument('--batchname', dest='BatchName', default="")
    parser.add_argument('--skimming_mode', action='store_true', default=False, help="Enable this option when anlyzer is skimmer.")
    parser.add_argument('--no_exec', action='store_true', default=False, help="only produce working area, not submitting to the condor pool")
    
    #Note: this option will change the behavior of the script. output directory will be changed to Your GV0, hadd will be disabled, and will create the info json of skimmed tree   
    return parser

def getMasterDirectoryName(timeStamp, Analyzer, Userflags):
    MasterDirectoryName = f"{timeStamp}_{Analyzer}"
    if len(Userflags) > 0:
        for flag in Userflags:
            MasterDirectoryName += f"_{flag}"
    abs_MasterDirectoryName = os.path.join(SKNANO_RUNLOG,MasterDirectoryName)
    print(f"Creating Master Working Directory: {abs_MasterDirectoryName}")
    #print(f"Copy library files to {MasterDirectoryName.split('/')}")
    #os.makedirs(abs_MasterDirectoryName)
    shutil.copytree(SKNANO_INSTALLDIR, abs_MasterDirectoryName) 
    #os.system(f"cp -r {SKNANO_INSTALLDIR} {abs_MasterDirectoryName}")
    #print("...Done")
    return MasterDirectoryName, abs_MasterDirectoryName

def getInputSampleList(inputArguments, eras):
    #if string
    if inputArguments.endswith(".txt"):
        with open(inputArguments, 'r') as f:
            InputSampleList = f.readlines()
    else:
        InputSampleList = inputArguments.split(",")
    InputSampleList = [x.replace("\n","").replace(" ","") for x in InputSampleList if not x.startswith("#") or not x.startswith("\n")]
    return InputSampleList
            
def jobFileDivider(files,ngroup):
    if ngroup > 0:
        if ngroup > len(files):
            ngroup = len(files)
        filegroups = [[] for i in range(ngroup)]
        for i, file in enumerate(files):
            filegroups[i%ngroup].append(file)
    elif ngroup < 0:
        #if ngroup is negative, it will be number of files per job
        nfile = abs(ngroup)
        if nfile > len(files):
            nfile = len(files)
        if nfile == 0:
            return []
        filegroups = [files[i:i + nfile] for i in range(0, len(files), nfile)]
        
    else:
        print('\033[91m'+"ERROR: ngroup should be positive or negative integer"+'\033[0m')
    return filegroups

def create_cpp_job(i, out_base, sample, argparse_module, sampleInfo, isMC, era, userflags, samplePaths, reduction, working_dir):
    """
    Worker function to create a single C++ job file.
    """
    output = out_base.replace('.root', f'_{i+1}.root')
    
    cpp_lines = [
        "#include <algorithm>",
        f"void job{i+1}() {{",
        f"    {argparse_module.Analyzer} module;",
        '    module.SetTreeName("Events");',
        '    module.LogEvery = 1000;',
    ]

    if isMC:
        cpp_lines += [
            '    module.IsDATA = false;',
            f'    module.MCSample = "{sample}";',
            f'    module.xsec = {sampleInfo["xsec"]};',
            f'    module.sumW = {sampleInfo["sumW"]};',
            f'    module.sumSign = {sampleInfo["sumsign"]};',
        ]
    else:
        cpp_lines += [
            '    module.IsDATA = true;',
            f'    module.DataStream = "{sample.split("_")[0]}";',
        ]

    cpp_lines += [
        f'    module.SetEra("{era}");',
    ]

    if userflags:
        cpp_lines.append("    module.Userflags = {")
        cpp_lines += [f'        "{flag}",' for flag in userflags]
        cpp_lines.append("    };")

    for path in samplePaths[i]:
        cpp_lines.append(f'    module.AddFile("{path}");')

    cpp_lines += [
        f'    module.SetOutfilePath("{output}");',
    ]

    if reduction > 1:
        cpp_lines.append(
            f'    module.MaxEvent = std::max(1, static_cast<int>(module.fChain->GetEntries()/{int(reduction)}));'
        )

    cpp_lines += [
        "    module.Init();",
        "    module.initializeAnalyzer();",
        "    module.Loop();",
        "    module.WriteHist();",
        "}",
    ]

    cpp_content = "\n".join(cpp_lines)
    job_filename = os.path.join(working_dir, f"job{i+1}.cpp")
    
    try:
        with open(job_filename, 'w') as f:
            f.write(cpp_content)
    except Exception as e:
        print(f"Error writing {job_filename}: {e}")

 
def pythonJobProducer(era, sample, argparse, masterJobDirectory, userflags, isample, totsamples):
    isMC, period = isMCandGetPeriod(sample)
    AnalyzerName = argparse.Analyzer

        
    working_dir = os.path.join(masterJobDirectory,era,sample)
        
    os.makedirs(working_dir)
    if SKIMMING_MODE:
        out_base, suffix = getSkimmingOutBaseAndSuffix(era, sample, AnalyzerName)
        if not os.path.exists(os.path.dirname(out_base)):
            os.makedirs(os.path.dirname(out_base))
    else:
        os.makedirs(os.path.join(working_dir,"output"))
        out_base = os.path.join(working_dir,"output","hists.root")
    njobs = argparse.NJobs
    reduction = argparse.Reduction
    
    if sample.startswith("Skim_"):
        SkimInfo = skimInfoJsons[era][sample if isMC else re.sub(f"_{period}$", "", sample)]
        sampleInfo = sampleInfoJsons[era][SkimInfo['PD']]
        samplePaths = json.load(open(os.path.join(SKNANO_DATA,era,'Sample','Skim',sample+'.json')))['path']
        sample = SkimInfo['PD']
    else:
        sampleInfo = sampleInfoJsons[era][sample if isMC else re.sub(f"_{period}$", "", sample)]
        samplePaths = json.load(open(os.path.join(SKNANO_DATA,era,'Sample','ForSNU',sample+'.json')))['path']
        
    samplePaths = jobFileDivider(samplePaths, njobs)

    
    totalNumberOfJobs = len(samplePaths)

    for i in tqdm(range(totalNumberOfJobs), position=1, leave=False, desc=f"Creating Jobs for {sample}, ({isample}/{totsamples})", smoothing=1.):
        output = out_base.replace('.root',f'_{i+1}.root')

        # Read the template file
        template_path = os.path.join(SKNANO_HOME, "templates", "job.cc")
        with open(template_path, 'r') as f:
            job_content = f.read()
            
        # Replace template variables
        job_content = job_content.replace("[jobname]", f"job_{i+1}")
        job_content = job_content.replace("[analyzer]", argparse.Analyzer)
        job_content = job_content.replace("[sample]", sample)
        job_content = job_content.replace("[era]", era)
        
        if isMC:
            job_content = job_content.replace("[xsec]", str(sampleInfo["xsec"]))
            job_content = job_content.replace("[sumW]", str(sampleInfo["sumW"]))
            job_content = job_content.replace("[sumSign]", str(sampleInfo["sumsign"]))
        else:
            # For data, remove MC-specific lines
            job_content = job_content.replace("module.IsDATA = false;", "module.IsDATA = true;")
            job_content = job_content.replace('module.MCSample = "[sample]";', f'module.DataStream = "{sample.split("_")[0]}";')
            job_content = job_content.replace("module.xsec = [xsec];", "")
            job_content = job_content.replace("module.sumW = [sumW];", "") 
            job_content = job_content.replace("module.sumSign = [sumSign];", "")

        # Handle userflags
        if userflags:
            userflags_str = "module.Userflags = {\n"
            userflags_str += "".join([f'\t"{flag}",\n' for flag in userflags])
            userflags_str += "    };"
        else:
            userflags_str = ""
        job_content = job_content.replace("[USERFLAGS]", userflags_str)

        # Handle sample paths
        samplepaths_str = "\n".join([f'\tmodule.AddFile("{path}");' for path in samplePaths[i]])
        job_content = job_content.replace("[SAMPLEPATHS]", samplepaths_str)

        # Handle reduction/max events
        maxevent_str = f'\tmodule.MaxEvent = std::max(1, static_cast<int>(module.fChain->GetEntries()/{int(reduction)}));'
        job_content = job_content.replace("[MAXEVENT]", maxevent_str)

        # Set output path
        job_content = job_content.replace("[output]", output)
        job_filename = os.path.join(working_dir, f"job_{i+1}.cc")
        with open(job_filename, 'w') as f:
            f.write(job_content)
            
    return working_dir, totalNumberOfJobs
            
def makeMainAnalyzerJobs(working_dir,abs_MasterDirectoryName,totalNumberOfJobs, argparse):
    nmax = argparse.NMax
    memory = argparse.Memory
    batchname = argparse.BatchName
    userflags = getUserFlagsList(argparse.Userflags)
    ncpu = argparse.ncpu
    if batchname == "":
        batchname = argparse.Analyzer
        if len(userflags) > 0:
            for flag in userflags:
                batchname += f"_{flag}"
    libpath = os.environ['LD_LIBRARY_PATH']
    libpath = libpath.split(":")
    libpath = [x for x in libpath if x != SKNANO_LIB]
    libpath = [os.path.join(abs_MasterDirectoryName,'lib')]+libpath
    libpath = ":".join(libpath)
    if 'ROOT_INCLUDE_PATH' not in os.environ:
        inclpath = ""
    else:
        inclpath = os.environ['ROOT_INCLUDE_PATH']
    inclpath = inclpath.split(":")
    inclpath = [x for x in inclpath if SKNANO_INSTALLDIR.split("/")[-1] not in x]
    inclpath = inclpath + [os.path.join(abs_MasterDirectoryName,'include')]
    inclpath = inclpath[-1] if len(inclpath) == 2 else ":".join(inclpath)
    
    template_path = os.path.join(SKNANO_HOME, "templates", "run.sh")
    with open(template_path, 'r') as f:
        run_content = f.read()
    run_content = run_content.replace("[SKNANO_HOME]", SKNANO_HOME)
    run_content = run_content.replace("[SKNANO_DATA]", SKNANO_DATA)
    run_content = run_content.replace("[WORKDIR]", working_dir)
    run_content = run_content.replace("[SKNANO_RUNLOG_LIB]", os.path.join(abs_MasterDirectoryName, 'lib'))
    #run_content = run_content.replace("[LD_LIBRARY_PATH]", libpath)
    run_content = run_content.replace("[ROOT_INCLUDE_PATH]", inclpath)
    with open(os.path.join(working_dir,"run.sh"),'w') as f:
        f.write(run_content)

    #submit condor jobs
    job_dict = main_job.copy()
    job_dict['JobBatchName'] = f"{batchname}_{working_dir.split('/')[-1]}_{working_dir.split('/')[-2]}"
    job_dict['executable'] = os.path.join(working_dir,"run.sh")
    job_dict['RequestMemory'] = f'ifthenelse(isUndefined(MemoryUsage),{memory},(MemoryUsage * 2))' # 2 times of memory usage
    job_dict['RequestCpus'] = ncpu
    job_dict['output'] = os.path.join(working_dir,"job_$(Process).out")
    job_dict['error'] = os.path.join(working_dir,"job_$(Process).err")
    job_dict['concurrency_limits'] = f"n{nmax}.{username}"
    
    return job_dict

def makeHaddJobs(working_dir,argparser,sample):
    AnalyzerName = argparser.Analyzer
    if len(userflags) > 0:
        for flag in userflags:
            AnalyzerName += f"_{flag}"
    era = working_dir.split('/')[-2]
    hadd_target = os.path.join(SKNANO_OUTPUT,AnalyzerName,era,sample+'.root')
    if not os.path.exists(os.path.dirname(hadd_target)):
        os.makedirs(os.path.dirname(hadd_target))

    template_path = os.path.join(SKNANO_HOME, "templates", "hadd.sh")
    with open(template_path, 'r') as f:
        hadd_content = f.read()
    hadd_content = hadd_content.replace("[WORKDIR]", working_dir)
    hadd_content = hadd_content.replace("[TARGET]", hadd_target)
    with open(os.path.join(working_dir,"hadd.sh"),'w') as f:
        f.write(hadd_content)
        
    job_dict = hadd_job.copy()
    job_dict['executable'] = os.path.join(working_dir,"hadd.sh")
    job_dict['JobBatchName'] = f"Hadd_{working_dir.split('/')[-1]}_{working_dir.split('/')[-2]}"
    job_dict['output'] = os.path.join(working_dir,"hadd.out")
    job_dict['error'] = os.path.join(working_dir,"hadd.err")

    return job_dict

def makeSkimPostProcsJobs(working_dir,sample, argparser,era):
    AnalyzerName = argparser.Analyzer
    isMC, period = isMCandGetPeriod(sample)
    out_base, suffix = getSkimmingOutBaseAndSuffix(era, sample, AnalyzerName) 
    out_base = os.path.dirname(out_base)
    if isMC:
        with open(os.path.join(working_dir,"postproc.sh"),'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(f"mv {out_base} {os.path.join(os.path.dirname(out_base),out_base.split('/')[-1].replace('Temp_',''))}\n")
            f.writelines(f"cd $SKNANO_PYTHON\n")
            f.writelines(f"python3 sampleManager.py --era {era} --makeSkimTreeInfo --skimTreeFolder {os.path.dirname(out_base)} --skimTreeSuffix {AnalyzerName.replace('Skim_','')} --skimTreeOrigPD {sample}\n")
    else:
        target_dir = os.path.join(os.path.dirname(os.path.dirname(out_base)),os.path.dirname(out_base).split("/")[-1].replace('Temp_',''))
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with open(os.path.join(working_dir,"postproc.sh"),'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(f"mv {out_base} {target_dir}\n")
            f.writelines(f"cd $SKNANO_PYTHON\n")
            f.writelines(f"python3 sampleManager.py --era {era} --makeSkimTreeInfo --skimTreeFolder {os.path.dirname(target_dir)} --skimTreeSuffix {AnalyzerName.replace('Skim_','')} --skimTreeOrigPD {sample}\n")
            f.writelines(f"""if [ -z "$(ls -A {os.path.dirname(out_base)})" ]; then\n""")
            f.writelines(f"\trmdir {os.path.dirname(out_base)}\n") 
            f.writelines(f"fi")
     
    job_dict = {}
    job_dict['executable'] = os.path.join(working_dir,"postproc.sh")
    job_dict['JobBatchName'] = f"PostProc_{working_dir.split('/')[-1]}_{working_dir.split('/')[-2]}"
    job_dict['universe'] = "vanilla"
    job_dict['getenv'] = "True"
    job_dict['RequestMemory'] = 1024
    job_dict['output'] = os.path.join(working_dir,"postproc.out")
    job_dict['error'] = os.path.join(working_dir,"postproc.err")
    job_dict['should_transfer_files'] = "YES"
    job_dict['when_to_transfer_output'] = "ON_EXIT"

    return job_dict
    


def getEachAnalyzerToPostDag(kwarg):
    analyzer_sub_dict = kwarg['analyzer_sub_dict']
    hadd_sub_dict = kwarg['hadd_sub_dict']
    totalNumberOfJobs = kwarg['totalNumberofJobs']
    batchname = kwarg['batchname']
    
    if totalNumberOfJobs == 0:
        return
    
    analyzer_layer = {
        'name' : f"Analyzer_{batchname}",
        'submit_description': htcondor.Submit(analyzer_sub_dict),
        'vars' : [{"Process":str(i)} for i in range(1,totalNumberOfJobs+1)]
    }
    
    hadd_layer = {
        'name' :  f"Postproc_{batchname}" if SKIMMING_MODE else f"Hadd_{batchname}",
        'submit_description' : htcondor.Submit(hadd_sub_dict)
    }
    
    return (analyzer_layer,hadd_layer)
        
def getFinalDag(hadd_layer_dicts,skim_postproc_layers,master_dir,argparser):
    batchname = argparser.BatchName
    userflags = getUserFlagsList(args.Userflags)
    if batchname == "":
        batchname = argparser.Analyzer
        batchname += "_".join(userflags) if userflags else ""
                
    dag_dir = os.path.join(master_dir,"dags")
    os.makedirs(dag_dir)

    job_dict = final_job.copy()
    job_dict['executable'] = os.path.join(dag_dir,"final.sh")
    job_dict['output'] = os.path.join(dag_dir,"final.out")
    job_dict['error'] = os.path.join(dag_dir,"final.err")

    final_content = os.path.join(SKNANO_HOME, "templates", "final.sh")  
    with open(final_content, 'r') as f:
        final_content = f.read()
    final_content = final_content.replace("[DAGDIR]", dag_dir)
    final_content = final_content.replace("[SKNANO_PYTHON]", os.environ['SKNANO_PYTHON'])
    final_content = final_content.replace("[TOKEN]", TOKEN)
    final_content = final_content.replace("[CHATID]", chat_id)
    final_content = final_content.replace("[MASTERDIR]", master_dir)
    with open(os.path.join(dag_dir,"final.sh"),'w') as f:
        f.write(final_content)
    
    dag = dags.DAG()
    finaldag = dags.DAG()
    hadd_layers = []
    postproc_layers = []
    if SKIMMING_MODE:
        for layer_dict in skim_postproc_layers:
            if layer_dict is None:
                continue
            analyzer_dict, postproc_dict = layer_dict
            if analyzer_dict is None or postproc_dict is None:
                continue
            analyzer_layer = dag.layer(
                name = analyzer_dict['name'],
                submit_description = analyzer_dict['submit_description'],
                vars = analyzer_dict['vars']
            )
            postproc_layer = analyzer_layer.child_layer(
                name = postproc_dict['name'],
                submit_description = postproc_dict['submit_description']
            )
            postproc_layers.append(postproc_layer)
        #if Skimming, each postproc layer will not run concurrently b/c io
        for i, layer in enumerate(postproc_layers):
            if i == 0:
                continue
            postproc_layers[i-1].add_parents(postproc_layers[i])

    else:
        for layer_dict in hadd_layer_dicts:
            if layer_dict is None:
                continue
            analyzer_dict, hadd_dict = layer_dict
            if analyzer_dict is None or hadd_dict is None:
                continue
            analyzer_layer = dag.layer(
                name = analyzer_dict['name'],
                submit_description = analyzer_dict['submit_description'],
                vars = analyzer_dict['vars']
            )
            hadd_layer = analyzer_layer.child_layer(
                name = hadd_dict['name'],
                submit_description = hadd_dict['submit_description']
            )
            hadd_layers.append(hadd_layer)
    


    
    dag_file = dags.write_dag(dag,dag_dir)
    subdag = finaldag.subdag(name='main_Dag',dag_file=str(dag_file))
    finalNode = finaldag.final(name=job_dict['JobBatchName'],submit_description=job_dict)
    with open(os.path.join(dag_dir,job_dict['JobBatchName']+'.sub'),'w') as f:
        for key, value in job_dict.items():
            f.writelines(f"{key} = {value}\n")
        f.writelines(f"queue\n")

    print(dag.describe())
    print(finaldag.describe())
    finalDag_file = dags.write_dag(finaldag,dag_dir,'finaldag.dag')

    if not args.no_exec:
        os.chdir(dag_dir)
        finalDag_submit = htcondor.Submit.from_dag(str(finalDag_file),{"force":1,"include_env":','.join(list(os.environ.keys())),"batch-name":batchname}) 
        cluster_id = htcondor.Schedd().submit(finalDag_submit).cluster()
        print(f"DAGMan job cluster is {cluster_id}")

    
if __name__ == '__main__':
    parser = setParser()
    args = parser.parse_args()
    eras = getEraList(args.Era, args.Run)
    SKIMMING_MODE = args.skimming_mode
    if args.Analyzer.startswith("Skim_") and not SKIMMING_MODE:
        print('\033[93m'+'''It seems like you want to skim the samples. If so, you need to enable skimming mode by passing the 
--skimming_mode arguments to write the output file to gv0 and generate the *.json file of the skimmed sample.\n'''+'\033[0m')
        #ask Y/N to continue, print in yellow color
        while True:
            answer = input('\033[93m'+"Do you want to enable the skimming mode? (Y/N): "+'\033[0m')
            if answer == 'Y' or answer == 'y':
                SKIMMING_MODE = True
                break
            elif answer == 'N' or answer == 'n':
                SKIMMING_MODE = False
                break
            else:
                print('\033[93m'+"Please enter Y or N"+'\033[0m')

    
    userflags = getUserFlagsList(args.Userflags)
    timestamp, string_JobStartTime = getTimeStamp()
    _, abs_MasterDirectoryName= getMasterDirectoryName(timestamp, args.Analyzer, userflags)
    InputSamplelist = getInputSampleList(args.InputSample, eras)
    
    dag_list = []
    hadd_layers = []
    postproc_layers = []

    for era in eras:
        print(f"Working on {era}")
        InputSamplelist_era = makeSampleList(InputSamplelist, era)
        for isample, sample in enumerate(InputSamplelist_era):
            working_dir, totalNumberofJobs = pythonJobProducer(era, sample, args, abs_MasterDirectoryName, userflags, isample, len(InputSamplelist_era))
            if totalNumberofJobs == None:
                continue
            analyzer_sub_dict = makeMainAnalyzerJobs(working_dir,abs_MasterDirectoryName,totalNumberofJobs,args)
            if SKIMMING_MODE:
                postproc_sub_dict = makeSkimPostProcsJobs(working_dir,sample,args,era)
            else:
                hadd_sub_dict = makeHaddJobs(working_dir,args,sample)
            
            if SKIMMING_MODE:
                dag_list.append({'era':era,'sample':sample,'analyzer_sub_dict':analyzer_sub_dict,'hadd_sub_dict':postproc_sub_dict,'totalNumberofJobs':totalNumberofJobs,'working_dir':working_dir,'batchname':f"{args.Analyzer}_{era}_{sample}"})
            else:
                dag_list.append({'era':era,'sample':sample,'analyzer_sub_dict':analyzer_sub_dict,'hadd_sub_dict':hadd_sub_dict,'totalNumberofJobs':totalNumberofJobs,'working_dir':working_dir,'batchname':f"{args.Analyzer}_{era}_{sample}"})
            if dag_list is not None:
                if SKIMMING_MODE:
                    postproc_layers.append(getEachAnalyzerToPostDag(dag_list[-1]))
                else:
                    hadd_layers.append(getEachAnalyzerToPostDag(dag_list[-1]))
            
    getFinalDag(hadd_layers, postproc_layers,abs_MasterDirectoryName, args)
