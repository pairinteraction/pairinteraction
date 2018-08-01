FROM ubuntu:bionic
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git file \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    clang clang-tidy \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

COPY boost.patch /tmp/
RUN cd /usr/include/boost && patch -p1 < /tmp/boost.patch

ENV CXX=clang++
