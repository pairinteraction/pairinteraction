#!/bin/sh

set -e;

cd "${SOURCE_DIR}";
mkdir -p build;
cd build;

if [ "$1" = "osx" ]; then
    cmake -DWITH_DMG=On -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
elif [ "$1" = "gcov" ]; then
    cmake -DWITH_COVERAGE=On -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
elif [ "$1" = "tidy" ]; then
    cmake -DWITH_CLANG_TIDY=On -DWITH_GUI=Off -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
else
    cmake -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
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
