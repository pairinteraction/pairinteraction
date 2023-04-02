FROM quay.io/pypa/manylinux2014_x86_64

RUN yum install -y \
    boost169-devel \
    cereal-devel \
    cmake3 \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    fontconfig \
    gsl-devel \
    libuuid-devel \
    openblas-devel \
    swig3

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
