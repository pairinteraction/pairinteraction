FROM opensuse:leap
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper -n --gpg-auto-import-keys install \
    patterns-openSUSE-devel_C_C++ gcc6-c++ \
    git cmake \
    python3-qt5-devel \
    sqlite3 sqlite3-devel \
    boost_1_61-devel \
    gsl-devel \
    swig python3 python3-devel \
    python3-numpy python3-numpy-devel python3-scipy \
    dpkg rpm-build \
  && zypper clean \
  && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-6 100 \
  && update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++-6 100
