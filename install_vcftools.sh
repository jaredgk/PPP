#!/bin/sh
set -ex
wget https://github.com/vcftools/vcftools/tarball/master -O vcftools.tar.gz
tar -xzvf vcftools.tar.gz
cd vcftools && ./configure --prefix=/usr && make && sudo make install
