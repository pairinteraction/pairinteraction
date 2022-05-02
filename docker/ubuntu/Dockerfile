FROM ubuntu:bionic
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git file \
    pyqt5-dev-tools \
    libboost-filesystem-dev libboost-program-options-dev \
    libboost-serialization-dev libboost-system-dev libboost-test-dev \
    libgsl-dev \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev \
    python3-numpy python3-scipy \
    dpkg-dev rpm \
    curl ca-certificates \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add(["CxxWrap", "SparseArrays", "Test", "SQLite"])'
