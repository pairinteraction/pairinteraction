FROM opensuse/leap:15.2
RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-devel-C-C++-devel_C_C++ \
    git cmake \
    python3-qt5-devel \
    sqlite3 sqlite3-devel \
    libboost_filesystem1_66_0-devel libboost_program_options1_66_0-devel \
    libboost_serialization1_66_0-devel libboost_system1_66_0-devel libboost_test1_66_0-devel \
    eigen3-devel gsl-devel \
    swig python3 python3-devel \
    python3-numpy python3-numpy-devel python3-scipy \
    dpkg rpm-build \
    curl ca-certificates \
  && zypper clean

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz \
  && tar xvf julia-1.0.5-linux-x86_64.tar.gz -C /usr --strip-components=1 \
  && rm julia-1.0.5-linux-x86_64.tar.gz \
  && JULIA_DEPOT_PATH=/usr/share/julia julia -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="CxxWrap",version="0.8.2"), Pkg.PackageSpec(name="SparseArrays"), Pkg.PackageSpec(name="Test"), Pkg.PackageSpec(name="SQLite",version="0.8.2")])'
