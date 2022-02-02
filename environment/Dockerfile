# hash:sha256:00459ceda7ce610bb15ce923ba7bc0c870364d1f760f3ed2372c3115a6c682dc
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake=3.10.2-1ubuntu2.18.04.2 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("nloptr")' \
    && Rscript -e 'remotes::install_version("survminer")'

COPY postInstall /
RUN /postInstall
