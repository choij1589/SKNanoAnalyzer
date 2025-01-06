#!/bin/sh
echo -e "\033[31m##################### WARNING ########################\033[0m"
echo -e "\033[31m####         THIS IS DEVELOPMENT VERSION          ####\033[0m"
echo -e "\033[31m######################################################\033[0m"
echo ""

# Set up environment
#export SKNANO_HOME=`pwd`
export SKNANO_HOME="/data9/Users/choij/SKNanoAnalyzer"
export SKNANO_RUNLOG="/gv0/Users/choij/SKNanoRunlog"
export SKNANO_OUTPUT="/gv0/Users/choij/SKNanoOutput"
echo "@@@@ Working Directory: $SKNANO_HOME"

# check configuration
CONFIG_FILE=$SKNANO_HOME/config/config.$USER
if [ -f "${CONFIG_FILE}" ]; then
    echo "@@@@ Reading configuration from $CONFIG_FILE"
    SYSTEM=$(grep '\[SYSTEM\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    PACKAGE=$(grep '\[PACKAGE\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export TOKEN_TELEGRAMBOT=$(grep '\[TOKEN_TELEGRAMBOT\]' "${CONFIG_FILE}" | cut -d' ' -f2)
    export USER_CHATID=$(grep '\[USER_CHATID\]' "${CONFIG_FILE}" | cut -d' ' -f2)

else
    echo "@@@@ Configuration file $CONFIG_FILE not found"
    echo "@@@@ Please create a configuration file in config/ with your username"
fi
echo "@@@@ System:  $SYSTEM"
echo "@@@@ Package: $PACKAGE"

# root configuration
# no cvmfs related configuration for conda
if [ $PACKAGE = "conda" ]; then
    echo "@@@@ Primary environment using conda"
    source /data9/Users/choij/miniconda3/bin/activate
    conda activate nano
elif [ $PACKAGE = "mamba" ]; then
    # set up mamba environment
    micronmamba activate Nano
elif [ $PACKAGE = "cvmfs" ]; then
    echo "@@@@ Primary environment using cvmfs"
    RELEASE="`cat /etc/redhat-release`"
    if [[ $RELEASE == *"7."* ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
    elif [[ $RELEASE == *"8."* ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos8-gcc12-opt/setup.sh
    elif [[ $RELEASE == *"9."* ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
    else
        echo "@@@@ Not running on redhat 7, 8, or 9"
        echo "@@@@ Consider using conda environment"
    fi
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
export SKNANO_BUILDDIR=$SKNANO_HOME/build
export SKNANO_INSTALLDIR=$SKNANO_HOME/install
export PATH=$SKNANO_PYTHON:$PATH
export PYTHONPATH=$PYTHONPATH:$SKNANO_PYTHON
export SKNANO_LIB=$SKNANO_INSTALLDIR/lib
export SKNANO_RUN3_NANOAODPATH="/gv0/DATA/SKNano/Run3NanoAODv12/"
export SKNANO_RUN2_NANOAODPATH="/gv0/DATA/SKNano/Run2NanoAODv9/"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_LIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_LIB

# setting LHAPDFs
if [ $PACKAGE = "conda" ]; then
    if [[ ! -d "external/lhapdf" ]]; then
        echo "@@@@ Installing LHAPDF for conda environment"
        ./scripts/install_lhapdf.sh
    fi
    export PATH=$PATH:$SKNANO_HOME/external/lhapdf/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/lib
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/lib
    export LHAPDF_DATA_PATH=$SKNANO_HOME/external/lhapdf/data
elif [ $PACKAGE = "cvmfs" ]; then
    echo "@@@@ configuring LHAPDF from cvmfs"
else
    echo "@@@@ LHAPDF not found"
fi
export LHAPDF_INCLUDE_DIR=`lhapdf-config --incdir`
export LHAPDF_LIB_DIR=`lhapdf-config --libdir`

echo "@@@@ LHAPDF include: $LHAPDF_INCLUDE_DIR"
echo "@@@@ LHAPDF lib: $LHAPDF_LIB_DIR"
echo "@@@@ reading data from $LHAPDF_DATA_PATH"

# setting up libtorch
if [[ ! -d "external/libtorch" ]]; then
    echo "@@@@ Installing LibTorch"
    ./scripts/install_libtorch.sh
fi
export LIBTORCH_INCLUDE_DIR=$SKNANO_HOME/external/libtorch/include
export LIBTORCH_LIB_DIR=$SKNANO_HOME/external/libtorch/lib
export LIBTORCH_INSTALL_DIR=$SKNANO_HOME/external/libtorch

# env for correctionlibs
if [ $PACKAGE = "conda" ]; then
    export CORRECTION_INCLUDE_DIR=`correction config --incdir`
    export CORRECTION_LIB_DIR=`correction config --libdir`
    export JSONPOG_INTEGRATION_PATH=$SKNANO_HOME/external/jsonpog-integration
elif [ $PACKAGE = "cvmfs" ]; then
    CORRECTION=`ll $(correction config --libdir | sed 's|/lib$||') | grep version.py | awk '{print $NF}' | sed 's|/lib/python3.9/site-packages/correctionlib/version.py||'`
    export CORRECTION_INCLUDE_DIR=${CORRECTION}/include
    export CORRECTION_LIB_DIR=${CORRECTION}/lib
    #export CORRECTION_INCLUDE_DIR=/cvmfs/sft.cern.ch/lcg/releases/correctionlib/2.2.2-c7cee/x86_64-el9-gcc13-opt/include
    #export CORRECTION_LIB_DIR=/cvmfs/sft.cern.ch/lcg/releases/correctionlib/2.2.2-c7cee/x86_64-el9-gcc13-opt/lib
    export JSONPOG_INTEGRATION_PATH=/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
else
    echo "@@@@ Correctionlib not found"
fi

#env for onnxruntime
if [ "$PACKAGE" = "conda" ]; then
    export ONNXRUNTIME_INCLUDE_DIR=${CONDA_PREFIX}/include/onnxruntime/core/session
    export ONNXRUNTIME_LIB_DIR=${CONDA_PREFIX}/lib
elif [ "$PACKAGE" = "cvmfs" ]; then
    export ONNXRUNTIME_INCLUDE_DIR=$(scram tool info onnxruntime | grep 'INCLUDE=' | awk -F= '{print $2}')
    export ONNXRUNTIME_LIB_DIR=$(scram tool info onnxruntime | grep 'LIBDIR=' | awk -F= '{print $2}')
else
    echo "@@@@ Onnxruntime not found"
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CORRECTION_LIB_DIR
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CORRECTION_LIB_DIR

echo "@@@@ Correction include: $CORRECTION_INCLUDE_DIR"
echo "@@@@ Correction lib: $CORRECTION_LIB_DIR"
