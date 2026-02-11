FROM debian:bookworm
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git git-lfs \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    locales \
    lcov \
    doxygen pandoc \
    ca-certificates \
    libgl1 \
    libgraphviz-dev \
    libxkbcommon-x11-0 freeglut3-dev \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*
