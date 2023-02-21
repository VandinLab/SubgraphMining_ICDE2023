#!/bin/sh
set -x;

# install nlopt locally
mkdir nlopt-2.7.1/build
cd nlopt-2.7.1/build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make
make install
cd ../..
mv nlopt-2.7.1/install/include/nlopt.h nlopt-2.7.1/install/include/mynlopt.h # to avoid conflicts

# compile code
cd src
if [ "$(uname)" == "Darwin" ]; then
  make all RPATH='-Xlinker -rpath ./../nlopt-2.7.1/install/lib'
else
  make all RPATH='-Wl,-rpath=./../nlopt-2.7.1/install/lib'
fi
cd ..
