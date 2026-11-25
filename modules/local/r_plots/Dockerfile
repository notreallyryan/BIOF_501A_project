FROM rocker/r-ver:4.4.1

# System libraries that ggtree and its dependencies need
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev procps \
    && rm -rf /var/lib/apt/lists/*

# R packages
RUN R -e "install.packages(c('BiocManager', 'ape', 'ggplot2', 'ggnewscale', 'optparse'))" \
 && R -e "BiocManager::install('ggtree', ask = FALSE, update = FALSE)"