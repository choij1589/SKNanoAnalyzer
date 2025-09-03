#!/bin/bash
HOSTNAME=`hostname`
echo -e "\033[32m##############################################################\033[0m"
echo -e "\033[32m                                                              \033[0m"
echo -e "\033[32m     ███████╗██╗  ██╗███╗   ██╗ █████╗ ███╗   ██╗ ██████╗     \033[0m"
echo -e "\033[32m     ██╔════╝██║ ██╔╝████╗  ██║██╔══██╗████╗  ██║██╔═══██╗    \033[0m"
echo -e "\033[32m     ███████╗█████╔╝ ██╔██╗ ██║███████║██╔██╗ ██║██║   ██║    \033[0m"
echo -e "\033[32m     ╚════██║██╔═██╗ ██║╚██╗██║██╔══██║██║╚██╗██║██║   ██║    \033[0m"
echo -e "\033[32m     ███████║██║  ██╗██║ ╚████║██║  ██║██║ ╚████║╚██████╔╝    \033[0m"
echo -e "\033[32m     ╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝     \033[0m"
echo -e "\033[32m                                                              \033[0m"
echo -e "\033[32m                       Version 1.0.0                          \033[0m"
echo -e "\033[32m##############################################################\033[0m"
echo ""

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
export SKNANO_HOME=`pwd`
export SKNANO_RUNLOG="/gv0/Users/$USER/SKNanoRunlog"
export SKNANO_OUTPUT="/data9/Users/$USER/SKNanoOutput"
echo "@@@@ Working Directory: $SKNANO_HOME"

CONFIG_FILE="$SKNANO_HOME/config/config.$USER"
# check configuration
if [ -f "${CONFIG_FILE}" ]; then
    echo -e "\033[32m@@@@ Reading configuration from $CONFIG_FILE\033[0m"
    PACKAGE=$(grep '\[PACKAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export TOKEN_TELEGRAMBOT=$(grep '\[TOKEN_TELEGRAMBOT\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export USER_CHATID=$(grep '\[USER_CHATID\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export SINGULARITY_IMAGE=$(grep '\[SINGULARITY_IMAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    # override SKNANO_* environment variables if they are set in the config file
    if grep -q '\[SKNANO_HOME\]' "${CONFIG_FILE}"; then
        export SKNANO_HOME=$(grep '\[SKNANO_HOME\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    fi
    if grep -q '\[SKNANO_RUNLOG\]' "${CONFIG_FILE}"; then
        export SKNANO_RUNLOG=$(grep '\[SKNANO_RUNLOG\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    fi
    if grep -q '\[SKNANO_OUTPUT\]' "${CONFIG_FILE}"; then
        export SKNANO_OUTPUT=$(grep '\[SKNANO_OUTPUT\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    fi
else
    echo -e "\033[31m@@@@ Configuration file $CONFIG_FILE not found\033[0m"
    echo -e "\033[3m@@@@ Please create a configuration file in config/ with your username\033[0m"
fi
echo "@@@@ System:  $SYSTEM"
echo "@@@@ Package: $PACKAGE"
echo "@@@@ Telegram Bot Token: $TOKEN_TELEGRAMBOT"
echo "@@@@ Telegram Chat ID:   $USER_CHATID"
echo "@@@@ Using singularity image: $SINGULARITY_IMAGE"
echo "@@@@ SKNano Home: $SKNANO_HOME"
echo "@@@@ SKNano Runlog: $SKNANO_RUNLOG"
echo "@@@@ SKNano Output: $SKNANO_OUTPUT"

# ROOT Package Settings
if [ $PACKAGE = "conda" ]; then
    echo -e "\033[32m@@@@ Primary environment using conda\033[0m"
    IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$IS_SINGULARITY" || -n "$GITHUB_ACTION" ]]; then
        # Building within Singularity image, will be used for batch jobs
        echo -e "\033[32m@@@@ Detected Singularity environment\033[0m"
        source /opt/conda/bin/activate
        conda activate torch
    else
        source ~/.conda-activate
        conda activate nano
    fi
elif [ $PACKAGE = "mamba" ]; then
    echo -e "\033[32m@@@@ Primary environment using mamba\033[0m"
    # Detect Singularity or GitHub Actions
    IS_SINGULARITY=$(env | grep -i "SINGULARITY_ENVIRONMENT")
    if [[ -n "$IS_SINGULARITY" || -n "$GITHUB_ACTION" ]]; then
        echo -e "\033[32m@@@@ Detected Singularity environment\033[0m"
    else
        # Only add micromamba to PATH if it's not already in PATH
        if ! command -v micromamba &> /dev/null; then
            export PATH="$HOME/micromamba/bin:${PATH}"
            export MAMBA_ROOT_PREFIX="$HOME/micromamba"
        fi
    fi

    # micromamba shell hook should only be run once per session
    if [[ -z "$__MAMBA_SETUP_DONE" ]]; then
        unalias mamba 2>/dev/null
        eval "$(micromamba shell hook -s zsh)"
        export __MAMBA_SETUP_DONE=1
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


export SKNANO_VERSION="Run3_v13_Run2_v9"
export SKNANO_DATA=$SKNANO_HOME/data/$SKNANO_VERSION
mkdir -p $SKNANO_DATA

export SKNANO_BIN=$SKNANO_HOME/bin
export SKNANO_PYTHON=$SKNANO_HOME/python
export SKNANO_BUILDDIR=$SKNANO_HOME/build/$SYSTEM
export SKNANO_INSTALLDIR=$SKNANO_HOME/install/$SYSTEM
export PATH=$SKNANO_PYTHON:$PATH
export PYTHONPATH=$PYTHONPATH:$SKNANO_PYTHON
export SKNANO_LIB=$SKNANO_INSTALLDIR/lib
export SKNANO_RUN3_NANOAODPATH="/gv0/DATA/SKNano/Run3NanoAODv13p1"
export SKNANO_RUN2_NANOAODPATH="/gv0/DATA/SKNano/Run2NanoAODv9p1"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_LIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_LIB

# setting LHAPDFs
if [[ ! -d "external/lhapdf/$SYSTEM" ]]; then
    echo -e "\033[32m@@@@ Installing LHAPDF for conda environment\033[0m"
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
    echo -e "\033[32m@@@@ Installing LibTorch\033[0m"
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

# ROCCOR
export ROCCOR_PATH=$SKNANO_HOME/external/RoccoR

# Check and copy RoccoR template files if missing
if [[ ! -f "$ROCCOR_PATH/CMakeLists.txt" ]] || [[ ! -f "$ROCCOR_PATH/RoccoR_LinkDef.hpp" ]]; then
    echo -e "\033[32m@@@@ Missing RoccoR build files, copying from templates...\033[0m"
    if [[ ! -f "$ROCCOR_PATH/CMakeLists.txt" ]]; then
        if [[ -f "$SKNANO_HOME/templates/RoccoR/CMakeLists.txt" ]]; then
            cp "$SKNANO_HOME/templates/RoccoR/CMakeLists.txt" "$ROCCOR_PATH/"
            echo "@@@@ Copied CMakeLists.txt to external/RoccoR/"
        else
            echo -e "\033[31m@@@@ Template CMakeLists.txt not found in templates/RoccoR/\033[0m"
        fi
    fi
    if [[ ! -f "$ROCCOR_PATH/RoccoR_LinkDef.hpp" ]]; then
        if [[ -f "$SKNANO_HOME/templates/RoccoR/RoccoR_LinkDef.hpp" ]]; then
            cp "$SKNANO_HOME/templates/RoccoR/RoccoR_LinkDef.hpp" "$ROCCOR_PATH/RoccoR_LinkDef.hpp"
            echo "@@@@ Copied RoccoR_LinkDef.hpp to external/RoccoR/RoccoR_LinkDef.hpp"
        else
            echo -e "\033[31m@@@@ Template RoccoR_LinkDef.hpp not found in templates/RoccoR/\033[0m"
        fi
    fi
fi

# JSONPOG integration auto-update
check_jsonpog_updates() {
    local auto_update=${1:-false}
    echo -e "\033[32m@@@@ Checking for updates in jsonpog-integration repository...\033[0m"
    export JSONPOG_REPO_PATH="$SKNANO_HOME/external/jsonpog-integration"

    if [ "$auto_update" = false ]; then
        echo -e "\033[32m@@@@ Auto-update is disabled. Skipping update check.\033[0m"
        return 0
    fi

    if [ ! -d "$JSONPOG_REPO_PATH" ]; then
        echo -e "\033[31m@@@@ JSONPOG Repository not found\033[0m"
        return 1
    fi
    
    cd "$JSONPOG_REPO_PATH"
    git fetch origin
    CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)

    if [ -z "$(git symbolic-ref -q HEAD)" ]; then
        echo -e "\033[32m@@@@ HEAD is detached. Switching back to the previous branch: $CURRENT_BRANCH\033[0m"
        git checkout "$CURRENT_BRANCH"
    fi

    # Check the current version of local and jsonpog version
    LOCAL_COMMIT_HASH=$(git rev-parse HEAD)
    LOCAL_COMMIT_DATE=$(git log -1 --format=%ci)
    echo -e "\033[32m@@@@ Current Local commit: $LOCAL_COMMIT_HASH\033[0m"
    echo -e "\033[32m@@@@ Current Local commit date: $LOCAL_COMMIT_DATE\033[0m"
    UPSTREAM_COMMIT_HASH=$(git ls-remote origin -h refs/heads/master | awk '{print $1}')
    UPSTREAM_COMMIT_DATE=$(git log -1 --format=%ci origin/master)
    echo -e "\033[32m@@@@ Latest JSONPOG (origin/master) commit: $UPSTREAM_COMMIT_HASH\033[0m"
    echo -e "\033[32m@@@@ Latest JSONPOG commit date: $UPSTREAM_COMMIT_DATE\033[0m"
    
    # Check if the local repository is behind the remote repository
    BEHIND=$(git rev-list --count origin/master..HEAD)

    if [ "$BEHIND" -gt 0 ]; then
        echo -e "\033[32m@@@@ Repository is $BEHIND commits behind origin/master.\033[0m"

        if [ "$auto_update" = true ]; then
            echo -e "\033[32m@@@@ Auto-update is enabled. Updating jsonpog-integration repository...\033[0m"
            git merge origin/master
            echo -e "\033[32m@@@@ Update completed!\033[0m"
        else
            echo -e "\033[32m@@@@ Auto-update is disabled. Skipping update.\033[0m"
        fi
    else
        echo -e "\033[32m@@@@ jsonpog-integration repository is already up-to-date.\033[0m"
    fi

    cd "$SKNANO_HOME"
}

# Call the function with auto_update set to false by default
# To enable auto-update, call with: check_jsonpog_updates true
check_jsonpog_updates false

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
