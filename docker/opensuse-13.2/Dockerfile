FROM opensuse:13.2
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
# jsoncpp is only in devel:tools
RUN zypper -n --gpg-auto-import-keys addrepo obs://build.opensuse.org/devel:tools devel:tools && \
  zypper -n --gpg-auto-import-keys refresh && \
  zypper -n --gpg-auto-import-keys install \
    patterns-openSUSE-devel_C_C++ \
    git cmake \
    eigen3-devel \
    python3-qt5-devel \
    lam lam-devel openmpi openmpi-devel \
    sqlite3 sqlite3-devel \
    boost-devel libboost_filesystem1_54_0 libboost_system1_54_0 \
    jsoncpp-devel \
    gsl-devel
