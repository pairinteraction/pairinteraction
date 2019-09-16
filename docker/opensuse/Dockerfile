FROM opensuse/leap:15.0
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper-n --gpg-auto-import-keys install --no-recommends \
    patterns-devel-C-C++-devel_C_C++ \
    git cmake \
    python3-qt5-devel \
    sqlite3 sqlite3-devel \
    libboost_*-devel \
    gsl-devel \
    swig python3 python3-devel \
    python3-numpy python3-numpy-devel python3-scipy \
    julia \
    dpkg rpm-build \
  && zypper clean
