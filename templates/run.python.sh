#!/bin/bash
export PATH=[MAMBA_BIN_PATH]:$PATH
export MAMBA_ROOT_PREFIX=[MAMBA_ROOT_PREFIX]
eval "$(micromamba shell hook -s bash)"
micromamba activate Nano

export SKNANO_HOME=[SKNANO_HOME]
export SKNANO_DATA=[SKNANO_DATA]
export SKNANO_BIN=$SKNANO_HOME/bin
export PYTHONPATH=$PYTHONPATH:$SKNANO_HOME/python

export PATH=$PATH:$SKNANO_HOME/external/lhapdf/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/redhat/lib
export LHAPDF_INCLUDE_DIR=$SKNANO_HOME/external/lhapdf/redhat/include
export LHAPDF_LIB_DIR=$SKNANO_HOME/external/lhapdf/redhat/lib
export CORRECTION_INCLUDE_DIR=`correction config --incdir`
export CORRECTION_LIB_DIR=`correction config --libdir`
export JSONPOG_REPO_PATH=$SKNANO_HOME/external/jsonpog-integration
export ROCCOR_PATH=$SKNANO_HOME/external/RoccoR
export ONNXRUNTIME_INCLUDE_DIR=${CONDA_PREFIX}/include/onnxruntime/core/session
export ONNXRUNTIME_LIB_DIR=${CONDA_PREFIX}/lib
export LIBTORCH_INCLUDE_DIR=$SKNANO_HOME/external/libtorch/include
export LIBTORCH_LIB_DIR=$SKNANO_HOME/external/libtorch/lib
export LIBTORCH_INSTALL_DIR=$SKNANO_HOME/external/libtorch

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CORRECTION_LIB_DIR:[SKNANO_RUNLOG_LIB]:$LIBTORCH_LIB_DIR

export ROOT_HIST=0
export ROOT_INCLUDE_PATH=[ROOT_INCLUDE_PATH]

cd [WORKDIR]

python job_$1.py
exit $?
