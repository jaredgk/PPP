#!/bin/sh
set -ex
wget https://faculty.washington.edu/browning/beagle/beagle.24Aug19.3e8.jar
mv beagle.24Aug19.3e8.jar beagle.jar
pip install -e ./
