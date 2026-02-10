#!/bin/bash
set -euo pipefail

INSTALL_DIR="$(pwd)/cvmh_install"

# set PATH_TO_AUTOMAKE to the actual path to your installed version of automake, something like "/usr/share/automake-1.16" 
PATH_TO_AUTOMAKE=""
if [ -z "$PATH_TO_AUTOMAKE" ]; then
    echo "Path to Automake has not been set! Set it first in 'generate_NetCDF.sh' and then run it again."
    exit 1
fi

python generate_cvmh_input.py

if [ -d ${INSTALL_DIR} ]; then
    echo "CVMH already installed, generating data..."
else 
    echo "Installing CVMH..."
    echo ${INSTALL_DIR} 
    mkdir ${INSTALL_DIR} 
    wget --no-check-certificate http://strike.scec.org/cvws/download/cvmh/cvmh-15.1.0.tar.gz
    tar zxvf cvmh-15.1.0.tar.gz
    cd cvmh-15.1.0

    cd aux-config
    
    ln -sf ${PATH_TO_AUTOMAKE}/depcomp
    ln -sf ${PATH_TO_AUTOMAKE}/install-sh
    ln -sf ${PATH_TO_AUTOMAKE}/missing
    cd ..

    echo ${PATH_TO_AUTOMAKE}/depcomp
    echo "Configuring..."
    ./configure --prefix=${INSTALL_DIR}
    make
    make install
    cd ..
    echo "Generating data..."
fi


${INSTALL_DIR}/bin/vx_lite -s -z dep -m ${INSTALL_DIR}/model < cvmh_in_off-fault > cvmh_out_off-fault
${INSTALL_DIR}/bin/vx_lite -s -z dep -m ${INSTALL_DIR}/model < cvmh_in_on-fault > cvmh_out_on-fault

python convert_cvmh_output.py
