FROM opensuse/leap:15.6 as builder

COPY cpp-httplib.spec /
RUN zypper -n --gpg-auto-import-keys refresh \
    && zypper -n --gpg-auto-import-keys install --no-recommends rpm-build \
    && rpmspec --query --buildrequires cpp-httplib.spec | xargs zypper -n install \
    && rpmbuild --undefine=_disable_source_fetch -bb --nocheck cpp-httplib.spec \
    && zypper clean

FROM opensuse/leap:15.6

WORKDIR /
COPY --from=builder /usr/src/packages/RPMS/x86_64/libcpp-httplib0_12-0.12.5-0.x86_64.rpm ./cpp-httlib.rpm
COPY --from=builder /usr/src/packages/RPMS/x86_64/cpp-httplib-devel-0.12.5-0.x86_64.rpm ./cpp-httlib-devel.rpm

RUN zypper -n --gpg-auto-import-keys refresh \
  && zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-devel-C-C++-devel_C_C++ gcc12-c++ \
    cmake tar gzip git git-lfs \
    doctest-devel spdlog-devel eigen3-devel fmt-devel openblas-devel \
    nlohmann_json-devel \
    doxygen pandoc \
    pkg-config \
    libglvnd \
    graphviz-devel \
  && zypper -n --gpg-auto-import-keys install --no-recommends --allow-unsigned-rpm \
    ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && rm ./cpp-httlib.rpm ./cpp-httlib-devel.rpm \
  && zypper clean
