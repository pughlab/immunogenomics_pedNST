# hash:sha256:10490210912887b7a572412f14761a5938388c243f6642bb44e442b0bbe797f3
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("nloptr")'

COPY postInstall /
RUN /postInstall
