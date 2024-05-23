FROM opensuse/leap:15.4 as build
COPY cpp-httplib.spec /
RUN zypper -n --gpg-auto-import-keys refresh \
    && zypper -n --gpg-auto-import-keys install --no-recommends rpm-build \
    && rpmspec --query --buildrequires cpp-httplib.spec | xargs zypper -n install \
    && rpmbuild --undefine=_disable_source_fetch -bb --nocheck cpp-httplib.spec

FROM opensuse/leap:15.4
COPY --from=build /usr/src/packages/RPMS/x86_64/libcpp-httplib0_12-0.12.5-0.x86_64.rpm ./cpp-httlib.rpm
COPY --from=build /usr/src/packages/RPMS/x86_64/cpp-httplib-devel-0.12.5-0.x86_64.rpm ./cpp-httlib-devel.rpm
RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-devel-C-C++-devel_C_C++ \
    git cmake \
    sqlite3 sqlite3-devel \
    libboost_headers1_66_0-devel \
    doctest-devel spdlog-devel cereal-devel eigen3-devel fmt-devel gsl-devel libuuid-devel openblas-devel \
    nlohmann_json-devel \
    swig python39 python39-devel libglvnd \
    dpkg rpm-build \
    curl \
    doxygen graphviz graphviz-devel pandoc \
  && zypper -n --gpg-auto-import-keys install --no-recommends --allow-unsigned-rpm \
    ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && rm ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && zypper clean

ADD https://bootstrap.pypa.io/pip/pip.pyz /pip.pyz
RUN python3.9 /pip.pyz --no-cache-dir install \
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
