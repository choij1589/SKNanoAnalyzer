#!/bin/bash
cd $SKNANO_HOME/external
echo -e "\033[33m@@@@ Downloading LibTorch\033[0m"
wget --no-verbose https://download.pytorch.org/libtorch/nightly/cpu/libtorch-shared-with-deps-latest.zip -O libtorch.zip
unzip -q libtorch.zip -d libtorch
cd libtorch
mv libtorch/* ./
rm -rf libtorch
cd ..
rm libtorch.zip
