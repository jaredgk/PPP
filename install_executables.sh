#!/bin/sh
set -ex
wget https://faculty.washington.edu/browning/beagle/beagle.24Aug19.3e8.jar
chmod +x beagle.24Aug19.3e8.jar
sudo mv beagle.24Aug19.3e8.jar /usr/local/bin/beagle.jar
#wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20191020.zip -O plink2.zip && unzip plink2.zip
#sudo mv plink2 /usr/local/bin
pip install -e ./
