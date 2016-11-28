FROM ubuntu:xenial
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install --no-install-recommends -y \
        cmake build-essential git \
        libeigen3-dev pyqt5-dev-tools \
        libboost-dev libboost-filesystem-dev libboost-system-dev \
        libgsl-dev \
        libopenmpi-dev openmpi-bin \
        libsqlite3-dev sqlite3 \
        dpkg rpm \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*
