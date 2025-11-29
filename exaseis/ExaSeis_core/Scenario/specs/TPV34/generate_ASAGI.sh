#!/bin/bash
set -euo pipefail

INSTALL_DIR="$(pwd)/cvmh_install"
mkdir ${INSTALL_DIR}

python generate_cvmh_input.py

wget --no-check-certificate http://strike.scec.org/cvws/download/cvmh/cvmh-15.1.0.tar.gz
tar zxvf cvmh-15.1.0.tar.gz
cd cvmh-15.1.0

cd aux-config
# replace "usr/share/automake-1.16" with actual path to your installed version of automake
ln -sf /usr/share/automake-1.16/depcomp
ln -sf /usr/share/automake-1.16/install-sh
ln -sf /usr/share/automake-1.16/missing
cd ..

./configure --prefix=${INSTALL_DIR}
make
make install
cd ..
${INSTALL_DIR}/bin/vx_lite -s -z dep -m ${INSTALL_DIR}/model < cvmh_in > cvmh_out

python convert_cvmh_output.py