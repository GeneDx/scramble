# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|

  config.vm.box = "bento/ubuntu-20.04"

  config.vm.hostname = "genedx"

  config.vm.provider "virtualbox" do |vb|
    vb.memory = 1024
    vb.cpus = 2
  end

  config.vm.provision "shell", inline: <<-SHELL
    apt-get -qq update
    apt-get install -q -y  \
      autoconf \
      autogen \
      build-essential \
      curl \
      libbz2-dev \
      libcurl4-openssl-dev \
      libhts-dev \
      liblzma-dev \
      libncurses5-dev \
      libnss-sss \
      libssl-dev \
      libxml2-dev \
      ncbi-blast+ \
      r-base \
      r-bioc-biostrings \
      r-bioc-rsamtools \
      r-cran-biocmanager \
      r-cran-devtools \
      r-cran-stringr \
      r-cran-optparse \
      zlib1g-dev
  SHELL
end
