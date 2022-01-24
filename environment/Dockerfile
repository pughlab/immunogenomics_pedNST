# hash:sha256:a94cf93a0e04774ff7517ad3c3207c8e06f18e317a298e31215877628578423c
FROM registry.codeocean.com/codeocean/r-base:3.6.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("ggbeeswarm", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gtable", "0.3.0")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")'
