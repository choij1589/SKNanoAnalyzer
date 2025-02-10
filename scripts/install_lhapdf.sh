#!/bin/bash
echo -e "\033[33mInstalling LHAPDF for $SYSTEM\033[0m"

LHAPDF_INSTALL_DIR=$SKNANO_HOME/external/lhapdf/$SYSTEM
LHAPDF_VERSION="6.5.5"
rm -rf $LHAPDF_INSTALL_DIR
mkdir -p $LHAPDF_INSTALL_DIR

# download and install lhapdf
cd $LHAPDF_INSTALL_DIR
wget --no-verbose https://lhapdf.hepforge.org/downloads/\?f\=LHAPDF-$LHAPDF_VERSION.tar.gz -O LHAPDF-$LHAPDF_VERSION.tar.gz
tar xf LHAPDF-$LHAPDF_VERSION.tar.gz
cd LHAPDF-$LHAPDF_VERSION

if [ $SKNANO_SYSTEM == "osx" ]; then
    ./configure --prefix=$LHAPDF_INSTALL_DIR CXX=clang++
else
    ./configure --prefix=$LHAPDF_INSTALL_DIR CXX=g++
fi
make -j8 && make install
cd $LHAPDF_INSTALL_DIR
rm -rf LHAPDF*

# download CMS central PDF sets
mkdir -p $SKNANO_HOME/external/lhapdf/data
cd $SKNANO_HOME/external/lhapdf/data
wget --no-verbose https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas.tar.gz
wget --no-verbose https://lhapdfsets.web.cern.ch/current/NNPDF31_nlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nnlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nlo_hessian_pdfas.tar.gz
cd $SKNANO_HOME
