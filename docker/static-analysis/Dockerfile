FROM ubuntu:jammy AS builder
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    python3 ca-certificates build-essential git cmake ninja-build libssl-dev \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

ARG DUCKDB_TAG=v1.1.3
RUN git clone --depth 1 --branch ${DUCKDB_TAG} --single-branch \
    https://github.com/duckdb/duckdb.git /tmp/duckdb
WORKDIR /tmp/duckdb
RUN GEN=ninja CMAKE_VARS_BUILD="-DBUILD_UNITTESTS=0 -DBUILD_SHELL=0" make release && \
    python3 scripts/amalgamation.py

FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir -p /opt/duckdb
COPY --from=builder /tmp/duckdb/build/release/src/libduckdb.so /opt/duckdb/libduckdb.so
COPY --from=builder /tmp/duckdb/src/amalgamation/duckdb.hpp /opt/duckdb/duckdb.hpp

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev liblapacke-dev libopenblas-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    locales \
    doxygen pandoc \
    clang clang-tidy clang-format \
    ca-certificates \
    libgl1 \
    libgraphviz-dev \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*
