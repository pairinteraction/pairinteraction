FROM debian:buster
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
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg rpm \
    lcov \
    doxygen graphviz ca-certificates \
    python3-sphinx python3-numpydoc python3-pip \
    jupyter-nbconvert pandoc python3-ipykernel \
    python3-matplotlib \
    intel-mkl-64bit-2019.1-053 \
  && curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])' \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/* \
  && pip3 install sphinx_rtd_theme
ENV MKLROOT="/opt/intel/compilers_and_libraries_2019/linux/mkl/"
