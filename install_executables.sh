#!/bin/sh
set -ex
wget https://github.com/vcftools/vcftools/tarball/master -O master_vcftools.tar.gz
tar -xzf master_vcftools.tar.gz
cd vcftools* && ./autogen.sh --silent && ./configure --prefix=/usr --silent && make --silent && sudo make install && cd ..
wget http://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar
mv beagle.21Jan17.6cc.jar beagle.jar
