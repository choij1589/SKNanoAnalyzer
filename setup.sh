#!/bin/bash
echo -e "\033[31m##################### WARNING ########################\033[0m"
echo -e "\033[31m####         THIS IS DEVELOPMENT VERSION          ####\033[0m"
echo -e "\033[31m######################################################\033[0m"
echo ""


# Set up environment
#export SKNANO_HOME=`pwd`
export SKNANO_HOME="/data9/Users/choij/Sync/workspace/SKNanoAnalyzer"
export SKNANO_RUNLOG="/gv0/Users/$USER/SKNanoRunlog"
export SKNANO_OUTPUT="/gv0/Users/$USER/SKNanoOutput"
echo -e "\033[33m@@@@ Working Directory: $SKNANO_HOME\033[0m"

# check os
if [[ "$(uname)" == "Darwin" ]]; then
    export SYSTEM="osx"
elif [[ -f "/etc/redhat-release" ]]; then
    export SYSTEM="redhat"
else
    echo "Unsupported OS"
    return 1
fi

# check configuration
CONFIG_FILE=$SKNANO_HOME/config/config.$USER
if [ -f "${CONFIG_FILE}" ]; then
    echo "@@@@ Reading configuration from $CONFIG_FILE"
    PACKAGE=$(grep '\[PACKAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export TOKEN_TELEGRAMBOT=$(grep '\[TOKEN_TELEGRAMBOT\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export USER_CHATID=$(grep '\[USER_CHATID\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export SINGULARITY_IMAGE=$(grep '\[SINGULARITY_IMAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
else
    echo "@@@@ Configuration file $CONFIG_FILE not found"
    echo "@@@@ Please create a configuration file in config/ with your username"
fi
echo -e "\033[33m@@@@ System:  $SYSTEM\033[0m"
echo -e "\033[33m@@@@ Package: $PACKAGE\033[0m"
echo -e "\033[33m@@@@ Telegram Bot Token: $TOKEN_TELEGRAMBOT\033[0m"
echo -e "\033[33m@@@@ Telegram Chat ID: $USER_CHATID\033[0m"
echo -e "\033[33m@@@@ Using singularity image: $SINGULARITY_IMAGE\033[0m"

# ROOT Package Settings
if [ $PACKAGE = "conda" ]; then
    echo "@@@@ Primary environment using conda"
    IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$IS_SINGULARITY" || -n "$GITHUB_ACTION" ]]; then
        # Building within Singularity image, will be used for batch jobs
        echo "@@@@ Detected Singularity environment"
        source /opt/conda/bin/activate
        conda activate torch
    else
        source ~/.conda-activate
        conda activate nano
    fi
elif [ $PACKAGE = "mamba" ]; then
    echo "@@@@ Primary environment using mamba"
    #IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$APPTAINER_NAME" || -n "$APPTAINER_NAME" ||  -n "$GITHUB_ACTION" ]]; then
        # Building within Singularity image, will be used for batch jobs
        echo "@@@@ Detected Singularity environment"
        export PATH="/opt/conda/bin:${PATH}"
        export MAMBA_ROOT_PREFIX="/opt/conda"
        eval "$(micromamba shell hook -s zsh)"
    else
        export PATH="$HOME/micromamba/bin:${PATH}"
        export MAMBA_ROOT_PREFIX="$HOME/micromamba"
        eval "$(micromamba shell hook -s zsh)"
    fi
    micromamba activate Nano
    # from this point on, we can follow conda version of setup
    PACKAGE="conda"
    alias conda="micromamba"
elif [ $PACKAGE = "cvmfs" ]; then
    echo -e "\033[31m@@@@ cvmfs is not supported anymore\033[0m"
    return 1
else
    echo "@@@@ Package not recognized"
    echo "@@@@ Please check configuration file in config/config.$USER"
fi
echo -e "\033[33m@@@@ ROOT path: $ROOTSYS\033[0m"


export SKNANO_VERSION="Run3_v12_Run2_v9"
export SKNANO_DATA=$SKNANO_HOME/data/$SKNANO_VERSION
mkdir -p $SKNANO_DATA

export SKNANO_BIN=$SKNANO_HOME/bin
export SKNANO_PYTHON=$SKNANO_HOME/python
export SKNANO_BUILDDIR=$SKNANO_HOME/build/$SYSTEM
export SKNANO_INSTALLDIR=$SKNANO_HOME/install/$SYSTEM
export PATH=$SKNANO_PYTHON:$PATH
export PYTHONPATH=$PYTHONPATH:$SKNANO_PYTHON
export SKNANO_LIB=$SKNANO_INSTALLDIR/lib
export SKNANO_RUN3_NANOAODPATH="/gv0/DATA/SKNano/Run3NanoAODv12/"
export SKNANO_RUN2_NANOAODPATH="/gv0/DATA/SKNano/Run2NanoAODv9/"
export ROOT_INCLUDE_DIRS=$ROOT_INCLUDE_DIRS:$SKNANO_INSTALLDIR/include

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_LIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_LIB

# setting LHAPDFs
if [[ ! -d "external/lhapdf/$SYSTEM" ]]; then
    echo "@@@@ Installing LHAPDF for conda environment"
    ./scripts/install_lhapdf.sh
    if [ $? -ne 0 ]; then
        echo -e "\033[31m@@@@ LHAPDF installation failed\033[0m"
        return 1
    fi
fi
export PATH=$PATH:$SKNANO_HOME/external/lhapdf/$SYSTEM/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/$SYSTEM/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/$SYSTEM/lib
export LHAPDF_DATA_PATH=$SKNANO_HOME/external/lhapdf/data
export LHAPDF_INCLUDE_DIR=$SKNANO_HOME/external/lhapdf/$SYSTEM/include
export LHAPDF_LIB_DIR=$SKNANO_HOME/external/lhapdf/$SYSTEM/lib
#export LHAPDF_INCLUDE_DIR=`lhapdf-config --incdir`
#export LHAPDF_LIB_DIR=`lhapdf-config --libdir`

echo -e "\033[33m@@@@ LHAPDF include: $LHAPDF_INCLUDE_DIR\033[0m"
echo -e "\033[33m@@@@ LHAPDF lib: $LHAPDF_LIB_DIR\033[0m"
echo -e "\033[33m@@@@ reading data from $LHAPDF_DATA_PATH\033[0m"

# setting up libtorch
if [[ ! -d "external/libtorch" ]]; then
    ./scripts/install_libtorch.sh
    if [ $? -ne 0 ]; then
        echo -e "\033[31m@@@@ LibTorch installation failed\033[0m"
        return 1
    fi
fi
export LIBTORCH_INCLUDE_DIR=$SKNANO_HOME/external/libtorch/include
export LIBTORCH_LIB_DIR=$SKNANO_HOME/external/libtorch/lib
export LIBTORCH_INSTALL_DIR=$SKNANO_HOME/external/libtorch

# env for correctionlib
CORRECTIONLIBS=$(conda list | grep "correctionlib")
if [ -z "$CORRECTIONLIBS" ]; then
    echo -e "\033[31m@@@@ correctionlib not found in conda environment\033[0m"
    echo -e "\033[31m@@@@ Please install correctionlib in conda environment\033[0m"
    return 1 
fi

export CORRECTION_INCLUDE_DIR=`correction config --incdir`
export CORRECTION_LIB_DIR=`correction config --libdir`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CORRECTION_LIB_DIR
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CORRECTION_LIB_DIR
echo -e "\033[33m@@@@ Correction include: $CORRECTION_INCLUDE_DIR"
echo -e "\033[33m@@@@ Correction lib: $CORRECTION_LIB_DIR"

# JSONPOG integration auto-update
export JSONPOG_INTEGRATION_PATH=$SKNANO_HOME/external/jsonpog-integration
echo -e "\033[33m@@@@ Checking for updates in jsonpog-integration repository...\033[0m"
cd "$JSONPOG_INTEGRATION_PATH"
git fetch origin 
LOCAL_HASH=$(git rev-parse HEAD) # local hash of latest
REMOTE_HASH=$(git rev-parse origin/master) # remote hash of latest
echo "@@@@ Local latest commit: $LOCAL_HASH"
echo "@@@@ Remote latest commit: $REMOTE_HASH"

# latest commit date and message of the remote branch
REMOTE_DATE=$(git log -1 --format=%ci origin/master) 
REMOTE_MESSAGE=$(git log -1 --format=%s origin/master) 
echo "@@@@ Remote latest update: $REMOTE_DATE - $REMOTE_MESSAGE"
cd -    
if [ "$LOCAL_HASH" != "$REMOTE_HASH" ]; then
    echo -e "\033[33m@@@@ Updating jsonpog-integration repository...\033[0m"
    git submodule update --remote external/jsonpog-integration
else
    echo -e "\033[33m@@@@ jsonpog-integration repository is already up-to-date.\033[0m"
fi

# env for onnxruntime
ONNXRUNTIME=$(conda list | grep "onnxruntime")
if [ -z "$ONNXRUNTIME" ]; then
    echo -e "\033[31m@@@@ onnxruntime not found in conda environment\033[0m"
    echo -e "\033[31m@@@@ Please install onnxruntime in conda environment\033[0m"
    return 1
fi
export ONNXRUNTIME_INCLUDE_DIR=${CONDA_PREFIX}/include/onnxruntime/core/session
export ONNXRUNTIME_LIB_DIR=${CONDA_PREFIX}/lib
echo -e "\033[33m@@@@ onnxruntime include: $ONNXRUNTIME_INCLUDE_DIR\033[0m"
echo -e "\033[33m@@@@ onnxruntime lib: $ONNXRUNTIME_LIB_DIR\033[0m"
