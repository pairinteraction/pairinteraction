FROM debian:stretch
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    build-essential cmake git \
    doxygen graphviz ca-certificates \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
