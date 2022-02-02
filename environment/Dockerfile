# hash:sha256:52f286af29117ce6e07e69feca528bb415cbd1a72141c95f8fdd6afe2f1681ee
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("nloptr")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap" \
    ))' # Original versions: 2.10.0

COPY postInstall /
RUN /postInstall
