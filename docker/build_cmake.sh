#!/bin/sh

set -e;

cd "${SOURCE_DIR}";
mkdir build;
cd build;

if [ "$1" = "win32" ]; then
    cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/mingw64.cmake ..;
elif [ "$1" = "osx" ]; then
    cmake -DWITH_DMG=On ..;
else
    cmake ..;
fi;

make;

if [ "$1" = "win32" ]; then
    make win32;
elif [ "$1" = "osx" ]; then
    make test;
    make package;
    make license;
else
    make test;
    make package;
fi;
