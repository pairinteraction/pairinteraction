FROM debian:buster
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    apt-transport-https ca-certificates curl gnupg2 \
  && curl -LO https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && echo 'deb http://apt.repos.intel.com/mkl all main' > /etc/apt/sources.list.d/intel-mkl.list \
  && apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git \
    pyqt5-dev-tools \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    julia \
    dpkg rpm \
    lcov \
    doxygen graphviz ca-certificates \
    python3-sphinx python3-numpydoc python3-pip \
    python3-nbconvert pandoc python3-ipykernel \
    python3-matplotlib \
    intel-mkl-64bit-2019.1-053 \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && pip3 install sphinx_rtd_theme
ENV MKLROOT="/opt/intel/compilers_and_libraries_2019/linux/mkl/"
