# hash:sha256:d4c33370c654c13e4fbabcdd0d088595cc3dcda42c90fc62a9289fb9e992388e
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake=3.10.2-1ubuntu2.18.04.2 \
        libnlopt-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.16")' \
    && Rscript -e 'remotes::install_version("circlize", "0.4.13")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("ggbeeswarm", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.1")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.3")' \
    && Rscript -e 'remotes::install_version("ggsignif", "0.6.3")' \
    && Rscript -e 'remotes::install_version("ggtext", "0.1.1")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gtable", "0.3.0")' \
    && Rscript -e 'remotes::install_version("nloptr")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-13")' \
    && Rscript -e 'remotes::install_version("survminer", "0.4.9")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap" \
    ))' # Original versions: 2.10.0

COPY postInstall /
RUN /postInstall
