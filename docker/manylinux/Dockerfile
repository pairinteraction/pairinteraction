FROM quay.io/pypa/manylinux2014_x86_64 as build

RUN yum install -y rpm-build
ADD spdlog.spec .
RUN yum-builddep -y spdlog.spec \
 && rpmbuild --undefine=_disable_source_fetch -bb spdlog.spec
ADD cpp-httplib.spec .
RUN yum-builddep -y cpp-httplib.spec \
 && rpmbuild --undefine=_disable_source_fetch -bb --without=tests cpp-httplib.spec

FROM quay.io/pypa/manylinux2014_x86_64

COPY --from=build /root/rpmbuild/RPMS/x86_64/spdlog-devel-1.7.0-1.el7.x86_64.rpm ./spdlog.rpm
COPY --from=build /root/rpmbuild/RPMS/x86_64/cpp-httplib-devel-0.12.5-1.el7.x86_64.rpm ./cpp-httplib.rpm

RUN yum install -y \
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
    swig3 \
    ./spdlog.rpm \
    ./cpp-httplib.rpm \
 && rm ./spdlog.rpm ./cpp-httplib.rpm

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
