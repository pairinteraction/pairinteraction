FROM debian:bullseye
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    apt-transport-https ca-certificates curl gnupg2 \
  && curl -LO https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
  && echo 'deb http://apt.repos.intel.com/mkl all main' > /etc/apt/sources.list.d/intel-mkl.list \
  && apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git \
    pyqt5-dev-tools \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    doctest-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev uuid-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    locales python3-psutil python3-pyqtgraph python3-pint \
    dpkg rpm \
    lcov \
    doxygen graphviz ca-certificates rsync \
    python3-sphinx python3-numpydoc python3-sphinx-rtd-theme \
    jupyter-nbconvert pandoc python3-ipykernel \
    python3-matplotlib \
    intel-mkl-64bit-2020.0-088 \
  && locale-gen --no-purge en_US.UTF-8 \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/*

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6

# https://github.com/boostorg/serialization/issues/219
COPY boost.patch /tmp/
RUN patch -d/ -p0 < /tmp/boost.patch

ENV MKLROOT="/opt/intel/compilers_and_libraries_2020/linux/mkl/"
