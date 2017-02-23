FROM opensuse:13.2
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper -n --gpg-auto-import-keys install \
    patterns-openSUSE-devel_C_C++ \
    git cmake \
    eigen3-devel \
    python3-qt5-devel \
    sqlite3 sqlite3-devel \
    boost-devel \
    gsl-devel \
    swig python3 python3-devel \
    dpkg rpm-build \
  && zypper clean
