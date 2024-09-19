FROM opensuse/leap:15.6 as build
COPY cpp-httplib.spec /
RUN zypper -n --gpg-auto-import-keys refresh \
    && zypper -n --gpg-auto-import-keys install --no-recommends rpm-build \
    && rpmspec --query --buildrequires cpp-httplib.spec | xargs zypper -n install \
    && rpmbuild --undefine=_disable_source_fetch -bb --nocheck cpp-httplib.spec

FROM opensuse/leap:15.6
COPY --from=build /usr/src/packages/RPMS/x86_64/libcpp-httplib0_12-0.12.5-0.x86_64.rpm ./cpp-httlib.rpm
COPY --from=build /usr/src/packages/RPMS/x86_64/cpp-httplib-devel-0.12.5-0.x86_64.rpm ./cpp-httlib-devel.rpm
COPY oneAPI.repo /etc/zypp/repos.d/
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
    pkg-config intel-oneapi-tbb-devel-2021.13 \
  && zypper -n --gpg-auto-import-keys install --no-recommends --allow-unsigned-rpm \
    ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && rm ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && zypper clean

ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest:$CMAKE_PREFIX_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin/:$PATH"
