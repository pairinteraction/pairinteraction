FROM ubuntu:focal
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git file \
    pyqt5-dev-tools \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    doctest-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev uuid-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    curl ca-certificates \
    doxygen graphviz \
    python3-sphinx python3-numpydoc python3-sphinx-rtd-theme \
    jupyter-nbconvert pandoc python3-ipykernel \
    python3-matplotlib \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/*

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6
