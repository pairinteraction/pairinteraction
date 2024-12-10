FROM debian:bookworm
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install --no-install-recommends -y \
    cmake build-essential git \
    doctest-dev libspdlog-dev libeigen3-dev libfmt-dev \
    libcpp-httplib-dev nlohmann-json3-dev \
    python3 python3-dev \
    locales \
    lcov \
    doxygen libgraphviz-dev pandoc \
    curl ca-certificates \
    libgl1 \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen \
  && apt-get clean -qq && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin/:$PATH"
