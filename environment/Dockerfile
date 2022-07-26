# hash:sha256:b9deab54b5b2582c87ff68f30d43ec04109c0a61d93e8482156701c899c675ea
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.4ubuntu1 \
        cmake=3.10.2-1ubuntu2.18.04.2 \
        curl \
        fftw3=3.3.7-1 \
        fftw3-dev=3.3.7-1 \
        ghostscript=9.26~dfsg+0-0ubuntu0.18.04.16 \
        libnlopt-dev=2.4.2+dfsg-4 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.18")' \
    && Rscript -e 'remotes::install_version("broom", "1.0.0")' \
    && Rscript -e 'remotes::install_version("circlize", "0.4.13")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("fftw", "1.0-6.1")' \
    && Rscript -e 'remotes::install_version("ggbeeswarm", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.1")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.3")' \
    && Rscript -e 'remotes::install_version("ggsignif", "0.6.3")' \
    && Rscript -e 'remotes::install_version("ggtext", "0.1.1")' \
    && Rscript -e 'remotes::install_version("grid")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gtable", "0.3.0")' \
    && Rscript -e 'remotes::install_version("nloptr", "2.0.3")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-13")' \
    && Rscript -e 'remotes::install_version("survminer", "0.4.9")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.2.0")'

RUN Rscript -e 'remotes::install_github( \
        "trinker/plotflow", \
        ref = "e21b80fcd604f6c2514b42a6f519453c85dc4fed")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap", \
        "qusage" \
    ))' # Original versions: 2.10.0 2.28.0

RUN mkdir /.code-server \
	&& cd /.code-server \
	&& curl -L "https://github.com/coder/code-server/releases/download/v4.0.2/code-server-4.0.2-linux-amd64.tar.gz" -o code-server.tar.gz \
	&& tar -xvf code-server.tar.gz \
	&& rm code-server.tar.gz \
	&& ln -s /.code-server/code-server-4.0.2-linux-amd64/bin/code-server  /usr/bin/code-server

COPY postInstall /
RUN /postInstall
