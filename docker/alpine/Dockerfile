FROM alpine:edge
RUN echo "http://dl-cdn.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories \
  && apk update \
  && apk add \
    bash cmake make g++ git file \
    eigen-dev py3-qt5 \
    boost-dev \
    gsl-dev \
    sqlite-dev sqlite \
    zeromq-dev \
    swig python3-dev py-numpy-dev \
    py3-numpy py3-scipy py3-zmq \
    dpkg rpm \
  && touch /usr/include/eigen3/signature_of_eigen3_matrix_library \
  && rm -rf /var/cache/apk/*
