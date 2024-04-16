FROM quay.io/pypa/manylinux2014_x86_64 as build

RUN yum install -y rpm-build
ADD spdlog.spec .
RUN yum-builddep -y spdlog.spec \
 && rpmbuild --undefine=_disable_source_fetch -bb spdlog.spec

FROM quay.io/pypa/manylinux2014_x86_64

COPY --from=build /root/rpmbuild/RPMS/x86_64/spdlog-devel-1.7.0-1.el7.x86_64.rpm ./spdlog.rpm

RUN yum install -y \
    boost169-devel \
    cereal-devel \
    doctest-devel \
    eigen3-devel \
    fmt-devel \
    fontconfig \
    gsl-devel \
    libuuid-devel \
    openblas-devel \
    swig3 \
    ./spdlog.rpm \
 && rm ./spdlog.rpm

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
