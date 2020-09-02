# source Image
FROM ubuntu:20.04

# set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt-get update and install global requirements
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y  \
      autoconf \
      autogen \
      build-essential \
      curl \
      libbz2-dev \
      libcurl4-openssl-dev \
      liblzma-dev \
      libncurses5-dev \
      libnss-sss \
      libssl-dev \
      libxml2-dev \
      ncbi-blast+ \
      r-base \
      zlib1g-dev

# apt-get clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install htslib
RUN curl -LO https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && \
    tar xfj htslib-1.10.2.tar.bz2 && \
    cd htslib-1.10.2 && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install

# install global r requirements
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(c('Biostrings'))"
RUN Rscript -e "library(devtools); install_github('mhahsler/rBLAST')"

# install scramble
COPY . /app
RUN cd /app/cluster_identifier/src && \
  make
RUN ln -s /app/cluster_identifier/src/build/cluster_identifier /usr/local/bin

# define default command
CMD ["Rscript"]
