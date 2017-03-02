#!/bin/sh
set -ex
wget https://github.com/vcftools/vcftools/tarball/master
wget https://github.com/vcftools/vcftools/tarball/master -O master_vcftools.tar.gz
tar -xzvf master_vcftools.tar.gz
cd vcftools* && ./configure --prefix=/usr && make && sudo make install
