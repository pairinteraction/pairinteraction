FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1

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
    graphviz-devel \
    && yum clean all

ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest:$CMAKE_PREFIX_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169

ENV PATH=/opt/python/cp39-cp39/bin/:$PATH
