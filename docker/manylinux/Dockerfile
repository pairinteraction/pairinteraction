FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1 AS builder

RUN yum install -y epel-release \
    && yum install -y \
    git cmake ninja-build curl perl-IPC-Cmd \
    && yum clean all

ENV OPENSSL_VERSION=3.0.16
WORKDIR /tmp/openssl3
RUN curl -fsSL "https://www.openssl.org/source/openssl-${OPENSSL_VERSION}.tar.gz" | tar -xz --strip-components=1
RUN ./config --prefix=/usr/local/openssl3 --openssldir=/usr/local/openssl3 --libdir=lib64 shared && \
    make -j"$(nproc)" && \
    make install_sw
ENV PATH="/usr/local/openssl3/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/openssl3/lib64:${LD_LIBRARY_PATH}"
ENV OPENSSL_ROOT_DIR="/usr/local/openssl3"

ARG DUCKDB_TAG=v1.1.3
RUN git clone --depth 1 --branch ${DUCKDB_TAG} --single-branch \
    https://github.com/duckdb/duckdb.git /tmp/duckdb
WORKDIR /tmp/duckdb
RUN GEN=ninja CMAKE_VARS_BUILD="-DBUILD_UNITTESTS=0 -DBUILD_SHELL=0" make release && \
    python3 scripts/amalgamation.py

FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1

COPY --from=builder /usr/local/openssl3 /usr/local/openssl3
ENV PATH="/usr/local/openssl3/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/openssl3/lib64:${LD_LIBRARY_PATH}"
ENV OPENSSL_ROOT_DIR="/usr/local/openssl3"

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
    graphviz-devel \
    libxkbcommon-x11 freeglut-devel \
    && yum clean all
