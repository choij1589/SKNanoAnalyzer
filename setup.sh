#!/bin/bash
echo -e "\033[31m##################### WARNING ########################\033[0m"
echo -e "\033[31m####         THIS IS DEVELOPMENT VERSION          ####\033[0m"
echo -e "\033[31m######################################################\033[0m"
echo ""
HOSTNAME=`hostname`
# check os
if [[ "$(uname)" == "Darwin" ]]; then
    export SYSTEM="osx"
elif [[ -f "/etc/redhat-release" ]]; then
    export SYSTEM="redhat"
else
    echo -e "\033[31mUnsupported OS\33[0m"
    return 1
fi

# Set up environment
if [[ $HOSTNAME == *"tamsa"* ]]; then
    export SKNANO_HOME="/data9/Users/choij/Sync/workspace/SKNanoAnalyzer"
    export SKNANO_RUNLOG="/gv0/Users/$USER/SKNanoRunlog"
    export SKNANO_OUTPUT="/data9/Users/choij/Sync/SKNanoOutput"
else
    export SKNANO_HOME=`pwd`
    export SKNANO_RUNLOG="/mnt/DUMP/$USER/SKNanoRunlog"
    export SKNANO_OUTPUT="/mnt/DUMP/$USER/SKNanoOutput"
fi
echo "@@@@ Working Directory: $SKNANO_HOME"

# check configuration
CONFIG_FILE=$SKNANO_HOME/config/config.$USER
if [ -f "${CONFIG_FILE}" ]; then
    echo -e "\033[33m@@@@ Reading configuration from $CONFIG_FILE\033[0m"
    PACKAGE=$(grep '\[PACKAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export TOKEN_TELEGRAMBOT=$(grep '\[TOKEN_TELEGRAMBOT\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export USER_CHATID=$(grep '\[USER_CHATID\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export SINGULARITY_IMAGE=$(grep '\[SINGULARITY_IMAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
else
    echo -e "\033[31m@@@@ Configuration file $CONFIG_FILE not found\033[0m"
    echo -e "\033[3m@@@@ Please create a configuration file in config/ with your username\033[0m"
fi
echo "@@@@ System:  $SYSTEM"
echo "@@@@ Package: $PACKAGE"
echo "@@@@ Telegram Bot Token: $TOKEN_TELEGRAMBOT"
echo "@@@@ Telegram Chat ID:   $USER_CHATID"
echo "@@@@ Using singularity image: $SINGULARITY_IMAGE"

# ROOT Package Settings
if [ $PACKAGE = "conda" ]; then
    echo -e "\033[33m@@@@ Primary environment using conda\033[0m"
    IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$IS_SINGULARITY" || -n "$GITHUB_ACTION" ]]; then
        # Building within Singularity image, will be used for batch jobs
        echo -e "\033[33m@@@@ Detected Singularity environment\033[0m"
        source /opt/conda/bin/activate
        conda activate torch
    else
        source ~/.conda-activate
        conda activate nano
    fi
elif [ $PACKAGE = "mamba" ]; then
    echo -e "\033[33m@@@@ Primary environment using mamba\033[0m"
    IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$IS_SINGULARITY" ||  -n "$GITHUB_ACTION" ]]; then
        # Building within Singularity image, will be used for batch jobs
        echo -e "\033[33m@@@@ Detected Singularity environment\033[0m"
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
echo "@@@@ ROOT path: $ROOTSYS"


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
export SKNANO_RUN3_NANOAODPATH="/gv0/Users/choij/SKNano"
export SKNANO_RUN2_NANOAODPATH="/gv0/Users/choij/SKNano"
#export SKNANO_RUN3_NANOAODPATH="/gv0/DATA/SKNano/Run3NanoAODv12/"
#export SKNANO_RUN2_NANOAODPATH="/gv0/DATA/SKNano/Run2NanoAODv9/"
#export ROOT_INCLUDE_DIRS=$ROOT_INCLUDE_DIRS:$SKNANO_INSTALLDIR/include

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_LIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_LIB

# setting LHAPDFs
if [[ ! -d "external/lhapdf/$SYSTEM" ]]; then
    echo -e "\033[33m@@@@ Installing LHAPDF for conda environment\033[0m"
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

echo "@@@@ LHAPDF include: $LHAPDF_INCLUDE_DIR"
echo "@@@@ LHAPDF lib: $LHAPDF_LIB_DIR"
echo "@@@@ reading data from $LHAPDF_DATA_PATH"

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
echo "@@@@ Correction include: $CORRECTION_INCLUDE_DIR"
echo "@@@@ Correction lib: $CORRECTION_LIB_DIR"

# JSONPOG integration auto-update
echo -e "\033[33m@@@@ Checking for updates in jsonpog-integration repository...\033[0m"
JSONPOG_REPO_PATH="$SKNANO_HOME/external/jsonpog-integration"

if [ ! -d "$JSONPOG_REPO_PATH" ]; then
    echo -e "\033[31m@@@@ JSONPOG Repository not found\033[0m"
else
    cd "$JSONPOG_REPO_PATH"
    git fetch origin
    CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)

    if [ -z "$(git symbolic-ref -q HEAD)" ]; then
        echo -e "\033[33m@@@@ HEAD is detached. Switching back to the previous branch: $CURRENT_BRANCH\033[0m"
        git checkout "$CURRENT_BRANCH"
    fi

    #Check the current version of local and jsonpog version
    LOCAL_COMMIT_HASH=$(git rev-parse HEAD)
    LOCAL_COMMIT_DATE=$(git log -1 --format=%ci)
    echo -e "\033[33m@@@@ Current Local commit: $LOCAL_COMMIT_HASH\033[0m"
    echo -e "\033[33m@@@@ Current Local commit date: $LOCAL_COMMIT_DATE\033[0m"
    UPSTREAM_COMMIT_HASH=$(git ls-remote origin -h refs/heads/master | awk '{print $1}')
    UPSTREAM_COMMIT_DATE=$(git log -1 --format=%ci origin/master)
    echo -e "\033[33m@@@@ Latest JSONPOG (origin/master) commit: $UPSTREAM_COMMIT_HASH\033[0m"
    echo -e "\033[33m@@@@ Latest JSONPOG commit date: $UPSTREAM_COMMIT_DATE\033[0m"
    
    # Check if the local repository is behind the remote repository
    BEHIND=$(git rev-list --count origin/master..HEAD)

    if [ "$BEHIND" -gt 0 ]; then
        echo -e "\033[33m@@@@ Repository is $BEHIND commits behind origin/master.\033[0m"

        # Check if the user wants to update
        while true; do
            read -p "Do you want to update jsonpog correction? (Y/N): " USER_INPUT
            case "$USER_INPUT" in
                [Yy]* ) 
                    echo -e "\033[33m@@@@ Updating jsonpog-integration repository...\033[0m"
                    git merge origin/master
                    echo -e "\033[33m@@@@ Update completed!\033[0m"
                    break  
                    ;;
                [Nn]* ) 
                    echo -e "\033[33m@@@@ Update skipped.\033[0m"
                    break  
                    ;;
                * ) 
                    echo -e "\033[31m@@@@ Invalid input. Please enter Y or N.\033[0m"
                    ;;
            esac
        done
    else
        echo -e "\033[33m@@@@ jsonpog-integration repository is already up-to-date.\033[0m"
    fi

    cd "$SKNANO_HOME"
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
echo "@@@@ onnxruntime include: $ONNXRUNTIME_INCLUDE_DIR"
echo "@@@@ onnxruntime lib: $ONNXRUNTIME_LIB_DIR"
