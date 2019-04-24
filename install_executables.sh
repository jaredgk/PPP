#!/bin/sh
set -ex
wget https://github.com/vcftools/vcftools/tarball/master -O master_vcftools.tar.gz
tar -xzf master_vcftools.tar.gz
cd vcftools* && ./autogen.sh --silent && ./configure --prefix=/usr --silent && make --silent && sudo make install && cd ..
wget http://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar
mv beagle.21Jan17.6cc.jar beagle.jar
wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 -O master_htslib.tar.bz2
tar -xjf master_htslib.tar.bz2
cd htslib* && ./configure --prefix=/usr --silent && make --silent && sudo make install && cd ..
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 -O master_bcftools.tar.bz2
tar -xjf master_bcftools.tar.bz2
cd bcftools* && ./configure --prefix=/usr --silent && make --silent && sudo make install && cd ../
pip install -e ./
