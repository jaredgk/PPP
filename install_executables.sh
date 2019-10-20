#!/bin/sh
set -ex
mkdir $TRAVIS_BUILD_DIR/bin
export PATH=$PATH:$TRAVIS_BUILD_DIR/bin/
wget https://faculty.washington.edu/browning/beagle/beagle.24Aug19.3e8.jar
mv beagle.24Aug19.3e8.jar $TRAVIS_BUILD_DIR/bin/beagle.jar
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20190826.zip -O plink2.zip && unzip plink2.zip
mv plink2 $TRAVIS_BUILD_DIR/bin
pip install -e ./
