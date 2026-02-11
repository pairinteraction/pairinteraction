FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git git-lfs \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev liblapacke-dev libopenblas-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    locales \
    doxygen pandoc \
    clang clang-tidy clang-format \
    ca-certificates \
    libgl1 \
    libgraphviz-dev \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*
