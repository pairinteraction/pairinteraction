FROM ubuntu:focal as build
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y \
    && apt-get install -y packaging-dev ubuntu-keyring devscripts equivs \
    && dget -ux http://archive.ubuntu.com/ubuntu/pool/universe/c/cpp-httplib/cpp-httplib_0.10.3+ds-1.dsc \
    && cd cpp-httplib-0.10.3+ds/ \
    && sed -i'' -e 's/debhelper-compat (= 13)/debhelper-compat (= 12)/g' debian/control \
    && yes | mk-build-deps --install --remove \
    && DEB_BUILD_OPTIONS=nocheck dpkg-buildpackage --build=binary --unsigned-changes \
    && cd - \
    && dget -ux http://archive.ubuntu.com/ubuntu/pool/universe/s/sphinxcontrib-mermaid/sphinxcontrib-mermaid_0.7.1-3.dsc \
    && cd sphinxcontrib-mermaid-0.7.1 \
    && sed -i'' -e 's/debhelper-compat (= 13)/debhelper-compat (= 12)/g' debian/control \
    && sed -i'' -e '/python3-myst-parser/d' debian/control \
    && yes | mk-build-deps --install --remove \
    && DEB_BUILD_OPTIONS=nocheck dpkg-buildpackage --build=binary --unsigned-changes

FROM ubuntu:focal
COPY --from=build /libcpp-httplib0_0.10.3+ds-1_amd64.deb ./libcpp-httplib.deb
COPY --from=build /libcpp-httplib-dev_0.10.3+ds-1_amd64.deb ./libcpp-httplib-dev.deb
COPY --from=build /python3-sphinxcontrib-mermaid_0.7.1-3_all.deb ./python3-sphinxcontrib-mermaid.deb
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git file \
    libboost-dev libboost-serialization-dev \
    doctest-dev libspdlog-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev liblapacke-dev libopenblas-dev uuid-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    clang-6.0 clang-tidy-6.0 \
    clang-format-6.0 python3-pip python3-setuptools black flake8 \
    curl ca-certificates netcat \
    doxygen graphviz \
    python3-sphinx python3-numpydoc python3-sphinx-rtd-theme \
    jupyter-nbconvert pandoc python3-ipykernel python3-nbsphinx \
    python3-matplotlib \
    ./libcpp-httplib.deb ./libcpp-httplib-dev.deb ./python3-sphinxcontrib-mermaid.deb \
  && apt-get clean -qq \
  && rm -rf /var/lib/apt/lists/* ./libcpp-httplib.deb ./libcpp-httplib-dev.deb ./python3-sphinxcontrib-mermaid.deb

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])' \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using CxxWrap, SparseArrays, Test, SQLite'

# https://github.com/JuliaLang/julia/issues/34276
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/julia/libstdc++.so.6

RUN ln -s /usr/bin/clang-6.0 /usr/local/bin/clang \
 && ln -s /usr/bin/clang++-6.0 /usr/local/bin/clang++ \
 && ln -s /usr/bin/clang-tidy-6.0 /usr/local/bin/clang-tidy \
 && ln -s /usr/bin/clang-format-6.0 /usr/local/bin/clang-format

ENV CXX=clang++
