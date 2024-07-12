FROM quay.io/pypa/manylinux_2_28_x86_64

RUN yum install -y epel-release \
    && yum install -y \
    boost169-devel \
    cereal-devel \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    json-devel \
    tbb-devel \
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

RUN cmake -DTBB_ROOT=/usr -DTBB_OS=Linux -DSAVE_TO=/usr/lib64/cmake/tbb -P /usr/lib64/cmake/tbb/tbb_config_generator.cmake

RUN curl -LsSf https://astral.sh/uv/install.sh | sh

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
