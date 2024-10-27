FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install --no-install-recommends -qq \
    cmake build-essential file \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev liblapacke-dev libopenblas-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    python3 python3-dev \
    doxygen libgraphviz-dev pandoc \
    curl ca-certificates \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin/:$PATH"
