FROM quay.io/pypa/manylinux2014_x86_64

RUN yum-config-manager --add-repo https://download.opensuse.org/repositories/science:/dlr/CentOS_7/science:dlr.repo \
 && yum install -y boost169-devel cereal-devel doctest-devel cmake3 eigen3-devel fmt-devel gsl-devel swig libuuid-devel fontconfig

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169
