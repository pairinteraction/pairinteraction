FROM quay.io/pypa/manylinux_2_28_x86_64

RUN yum install -y epel-release \
    && yum install -y \
    boost169-devel \
    cereal-devel \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    json-devel \
    fontconfig \
    gsl-devel \
    libuuid-devel \
    openblas-devel \
    swig \
    spdlog-devel \
    cpp-httplib-devel

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
