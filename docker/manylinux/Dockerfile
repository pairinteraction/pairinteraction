FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1 AS builder

RUN yum install -y epel-release \
    && yum install -y \
    git cmake ninja-build openssl-devel \
    && yum clean all

ARG DUCKDB_TAG=v1.1.3
RUN git clone --depth 1 --branch ${DUCKDB_TAG} --single-branch \
    https://github.com/duckdb/duckdb.git /tmp/duckdb
WORKDIR /tmp/duckdb
RUN GEN=ninja CMAKE_VARS_BUILD="-DBUILD_UNITTESTS=0 -DBUILD_SHELL=0" make release && \
    python3 scripts/amalgamation.py

FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1

RUN mkdir -p /opt/duckdb
COPY --from=builder /tmp/duckdb/build/release/src/libduckdb.so /opt/duckdb/libduckdb.so
COPY --from=builder /tmp/duckdb/src/amalgamation/duckdb.hpp /opt/duckdb/duckdb.hpp

RUN yum install -y epel-release \
    && yum install -y \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    json-devel \
    openblas-devel \
    spdlog-devel \
    cpp-httplib-devel \
    libglvnd \
    && yum clean all

ENV PATH=/opt/python/cp39-cp39/bin/:$PATH
