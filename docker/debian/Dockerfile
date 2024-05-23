FROM debian:bookworm
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    apt-transport-https ca-certificates curl \
  && curl -sSfLo /etc/apt/trusted.gpg.d/oneapi-archive-keyring.asc https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  && echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list \
  && apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git \
    libboost-dev \
    doctest-dev libspdlog-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev uuid-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev libgl1 \
    locales \
    dpkg rpm \
    lcov \
    doxygen graphviz libgraphviz-dev ca-certificates rsync pandoc \
    intel-oneapi-mkl-devel \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/*

ADD https://bootstrap.pypa.io/pip/pip.pyz /pip.pyz
RUN python3 /pip.pyz --no-cache-dir install --break-system-packages \
    scikit-build-core \
    numpy \
    scipy \
    matplotlib \
    PyQt5 \
    pyqtgraph \
    Pint \
    pydantic \
    pytest \
    Sphinx \
    numpydoc \
    sphinx-rtd-theme \
    sphinxcontrib-mermaid \
    autodoc_pydantic[erdantic] \
    nbsphinx \
    jupyter \
    nbconvert
RUN rm -f /pip.pyz

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])' \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using CxxWrap, SparseArrays, Test, SQLite'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6
