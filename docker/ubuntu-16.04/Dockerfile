FROM ubuntu:xenial
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git \
    libeigen3-dev pyqt5-dev-tools \
    libboost-all-dev \
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    libzmq3-dev \
    swig python3 python3-dev \
    dpkg rpm \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
