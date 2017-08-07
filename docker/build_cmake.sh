#!/bin/sh

set -e;

cd "${SOURCE_DIR}";
mkdir build;
cd build;

if [ "$1" = "osx" ]; then
    cmake -DWITH_DMG=On ..;
elif [ "$1" = "gcov" ]; then
    cmake -DWITH_COVERAGE=On ..;
else
    cmake ..;
fi;

make -j 2;

if [ "$1" = "osx" ]; then
    make check;
    make package;
    make license;
else
    make check;
    make package;
fi;
