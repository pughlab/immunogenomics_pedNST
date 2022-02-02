# hash:sha256:003be8b1b3fa8728387ec5808634041b4766c9763aba76f8b76885097b72fc9d
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("circlize", "0.4.13")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("ggbeeswarm", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.3")' \
    && Rscript -e 'remotes::install_version("ggsignif", "0.6.3")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gtable", "0.3.0")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-13")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap", \
        "survminer" \
    ))' # Original versions: 2.10.0 latest

COPY postInstall /
RUN /postInstall
