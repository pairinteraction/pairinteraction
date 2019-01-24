FROM ubuntu:bionic
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    wget ca-certificates gnupg2 \
  && wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list \
  && apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git file \
    pyqt5-dev-tools \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    intel-mkl-64bit-2019.1-053 \
  && apt-get purge wget ca-certificates gnupg2 \
  && apt-get autoremove --purge \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
ENV MKLROOT="/opt/intel/compilers_and_libraries_2019/linux/mkl/"
