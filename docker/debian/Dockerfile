FROM debian:bookworm
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y wget gnupg ca-certificates
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git \
    libboost-dev \
    doctest-dev libspdlog-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev uuid-dev \
    libcpp-httplib-dev nlohmann-json3-dev intel-oneapi-tbb-devel-2021.13 \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev libgl1 \
    locales \
    dpkg rpm \
    lcov \
    curl ca-certificates \
    doxygen graphviz libgraphviz-dev ca-certificates rsync pandoc \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/*

ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest:$CMAKE_PREFIX_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

ADD https://bootstrap.pypa.io/pip/pip.pyz /pip.pyz
RUN python3 /pip.pyz --no-cache-dir install --break-system-packages \
    scikit-build-core \
    numpy \
    scipy \
    matplotlib \
    mkl \
    mkl-devel \
    nanobind \
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

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin/:$PATH"

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])' \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using CxxWrap, SparseArrays, Test, SQLite'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6

# For compatiblity (to be removed)
RUN mkdir -p /opt/intel/oneapi/ && touch /opt/intel/oneapi/setvars.sh
