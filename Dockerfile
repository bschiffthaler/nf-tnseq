ARG R_VERSION=4.0.2

FROM r-base:${R_VERSION}

RUN apt update && \
    apt install -y libxml2-dev libcurl4-openssl-dev libssl-dev pandoc

RUN Rscript -e "install.packages(c('Rcpp', 'tidyverse', 'pander', 'stringdist', \
                                   'plotly', 'BiocManager', 'readxl', \
                                   'writexl', 'here', 'rmarkdown', 'reactable', \
                                   'corrplot', 'ggrepel'),\
                                 Ncpus = 6)"

RUN Rscript -e "BiocManager::install(c('Tnseq', 'Gviz', 'Rsubread', 'pheatmap'), Ncpus = 6)"

ENTRYPOINT ["Rscript"]