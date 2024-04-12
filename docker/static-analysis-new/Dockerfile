FROM ubuntu:focal
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git file \
    libboost-dev libboost-serialization-dev \
    doctest-dev libspdlog-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev liblapacke-dev libopenblas-dev uuid-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    clang clang-tidy \
    clang-format python3-pip python3-setuptools black flake8 \
    curl ca-certificates netcat \
    doxygen graphviz \
    python3-sphinx python3-numpydoc python3-sphinx-rtd-theme \
    jupyter-nbconvert pandoc python3-ipykernel python3-nbsphinx \
    python3-matplotlib \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && apt-get install -y debhelper dh-python \
 && mkdir -p /src && cd /src \
 && curl -LOOO http://archive.ubuntu.com/ubuntu/pool/universe/s/sphinxcontrib-mermaid/sphinxcontrib-mermaid_0.7.1{-3.dsc,.orig.tar.gz,-3.debian.tar.xz} \
 && dpkg-source -x sphinxcontrib-mermaid_0.7.1-3.dsc \
 && cd sphinxcontrib-mermaid-0.7.1 \
 && DEB_BUILD_OPTIONS=nocheck dpkg-buildpackage -b -d -uc -us \
 && cd - \
 && apt-get install ./python3-sphinxcontrib-mermaid_0.7.1-3_all.deb \
 && cd / && rm -rf /src \
 && apt-get purge --autoremove -y debhelper dh-python \
 && apt-get clean -qq \
 && rm -rf /var/lib/apt/lists/*

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])' \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using CxxWrap, SparseArrays, Test, SQLite'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6
