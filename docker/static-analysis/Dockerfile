FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y curl gnupg ca-certificates \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*
RUN curl -LsSf https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git file \
    libboost-dev libboost-serialization-dev \
    doctest-dev libspdlog-dev libcereal-dev libeigen3-dev libfmt-dev libgsl-dev liblapacke-dev libopenblas-dev uuid-dev \
    libcpp-httplib-dev nlohmann-json3-dev intel-oneapi-tbb-devel-2021.13 \
    libsqlite3-dev sqlite3 \
    swig python3 python3-dev libgl1 \
    locales \
    dpkg-dev rpm \
    ca-certificates netcat \
    doxygen graphviz libgraphviz-dev pandoc \
    clang clang-tidy clang-format \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest:$CMAKE_PREFIX_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin/:$PATH"
