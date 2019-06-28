#!/bin/sh
set -ex
wget http://cran.rstudio.com/src/base/R-3/R-3.4.1.tar.gz
tar xvf R-3.4.1.tar.gz
cd R-3.4.1 && ./configure --prefix=/usr/R --silent && make && sudo make install && cd ..