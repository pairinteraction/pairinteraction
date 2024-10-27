FROM quay.io/pypa/manylinux_2_28_x86_64:2024.09.16-1

RUN yum install -y epel-release \
    && yum install -y \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    json-devel \
    openblas-devel \
    spdlog-devel \
    cpp-httplib-devel \
    graphviz-devel \
    && yum clean all

ENV PATH=/opt/python/cp39-cp39/bin/:$PATH
