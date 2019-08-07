#!/bin/sh
set -ex
wget http://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar
mv beagle.21Jan17.6cc.jar beagle.jar
pip install -e ./
