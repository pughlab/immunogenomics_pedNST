# hash:sha256:7180d8b094836e1927a46802c0c308df531ca8a71d7c10c790f85ce5af482ca5
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("nloptr")'

COPY postInstall /
RUN /postInstall
