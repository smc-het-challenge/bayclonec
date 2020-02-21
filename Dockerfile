FROM r-base:3.5.0

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('mclust','fpc','DPpackage'), repos='http://cran.rstudio.com/')"


RUN apt-get install -y git

WORKDIR /opt

RUN git clone https://github.com/compgenome365/bayclonec.git
 
WORKDIR /opt/bayclonec
 
RUN g++ -o parseInputData_smc parseInputData_smc.cpp
 
ENTRYPOINT ["/opt/bayclonec/run_bc.sh"]
