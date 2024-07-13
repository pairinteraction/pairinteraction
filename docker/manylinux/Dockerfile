FROM quay.io/pypa/manylinux_2_28_x86_64

COPY oneAPI.repo /etc/yum.repos.d/oneAPI.repo

RUN yum install -y epel-release \
    && yum install -y \
    boost169-devel \
    cereal-devel \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    json-devel \
    intel-oneapi-tbb-devel-2021.13 \
    fontconfig \
    gsl-devel \
    libuuid-devel \
    openblas-devel \
    swig \
    spdlog-devel \
    cpp-httplib-devel \
    python3 \
    python3-devel \
    graphviz-devel

ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest:$CMAKE_PREFIX_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

RUN curl -LsSf https://astral.sh/uv/install.sh | sh

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
