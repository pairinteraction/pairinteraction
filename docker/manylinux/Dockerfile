FROM quay.io/pypa/manylinux2014_x86_64

RUN yum install -y boost169-devel cmake3 eigen3-devel gsl-devel swig3

ENV BOOST_INCLUDEDIR=/usr/include/boost169
ENV BOOST_LIBRARYDIR=/usr/lib64/boost169

# Use Python 3.8
ENV PATH="/opt/python/cp38-cp38/bin/:$PATH"
ENV PYTHON_INCLUDE_DIR="/opt/python/cp38-cp38/include/python3.8/"

# TODO: Install Python packages from requirements.txt
RUN pip3 install --only-binary=:all: numpy scipy pyqt5
