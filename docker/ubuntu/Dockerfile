FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential git git-lfs file \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev liblapacke-dev libopenblas-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    doxygen pandoc \
    ca-certificates \
    libgl1 \
    libgraphviz-dev \
    libxkbcommon-x11-0 freeglut3-dev \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*
