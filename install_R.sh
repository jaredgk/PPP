#!/bin/sh
set -ex
sudo apt-get update
sudo apt-get install gfortran
wget http://cran.rstudio.com/src/base/R-3/R-3.4.1.tar.gz
tar -xf R-3.4.1.tar.gz
cd R-3.4.1 && ./configure --prefix=/usr/R --enable-R-shlib --silent && make --silent && sudo make install && cd ..
export PATH=/usr/R/bin:$PATH
export R_HOME=/usr/R/lib64/R
whereis R

