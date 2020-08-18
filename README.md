Scramble
========

Repository Layout
-----------------
1. `cluster_identifier/` - this directory contains the part of the application responsible for identifying soft clipped
clusters. For how to build see the build section.
2. `cluster_analysis/` - contains code related to the analysis of soft-clipped clusters and interpretation.
3. `validation/` - sample bam and outputs for testing SCRAMble
4. `README` - Repository description/notes.

Build
-----

Install dependencies (Debian/Ubuntu):

    $ apt-get update
    $ apt-get install -y  \
        build-essential \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libnss-sss \
        libssl-dev \
        libhts-dev \
        r-base \
        zlib1g-dev

Install dependencies (Centos):

    $ yum install -y \
          epel-release && \
    $ yum install -y \
          autoconf \
          bzip2 \
          bzip2-devel \
          libcurl-devel \
          libcrypto-devel \
          openssl-devel \
          R \
          xz-devel \
          zlib-devel

Install R packages dependencies:

    $ Rscript -e "install.packages('optparse')"
    $ Rscript -e "install.packages('stringr')"
    $ Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biostrings')"

To build the cluster_identifier (estimated install time <5 minutes):

    $ cd cluster_identifier/src
    $ make

That should be it. It will create an executable named `build/cluster_identifier`.
 
Building requires `HTSlib` and a few other dev packages (installation instructions for Debian/Ubuntu/Centos above).
Please edit `cluster_identifier/src/Makefile` if `HTSlib` is not installed in the default location.

    /usr/local/include/htslib/*.h
    /usr/local/lib/libhts.a

Running
-------
SCRAMble runs as a two-step process. First `cluster_identifier` is used to generate soft-clipped read cluster consensus
sequences. Second, `SCRAMble-MEIs.R` analyzes the cluster file for likely MEIs. Running SCRAMble on the test bam in the validation directory should take <1 minute for each step.

To run SCRAMble cluster_identifier:

    $ /path/to/scramble/cluster_identifier/src/build/cluster_identifier \
        /path/to/install_dir/scramble/validation/test.bam > /path/to/output/clusters.txt

To run SCRAMble-MEIs and SCRAMble-dels(with default settings):

    $ Rscript --vanilla /path/to/scramble/cluster_analysis/bin/SCRAMble-MEIs.R \
        --out-name /path/to/output/out 	\
        --cluster-file /path/to/output/clusters.txt \
        --install-dir /path/to/scramble/cluster_analysis/bin \
        --mei-refs /path/to/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa \
        --blast-ref /path/to/scramble/cluster_analysis/resources/blast_refs/chr2.fa \
        --eval-meis \
        --eval-dels 
	
Running with Docker
-------------------
SCRAMble is also distributed with a `Dockerfile`. Running SCRAMble using `docker` (estimated install time <10 minutes):

    $ git clone https://github.com/GeneDx/scramble.git
    $ cd scramble
    $ docker build -t scramble:latest .
    $ docker run -it --rm scramble:latest bash
    # cluster_identifier \
        /app/validation/test.bam > clusters.txt
    # Rscript --vanilla /app/cluster_analysis/bin/SCRAMble-MEIs.R \
        --out-name ${PWD}/out \
        --cluster-file ${PWD}/clusters.txt \
        --install-dir /app/cluster_analysis/bin \
        --mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
        --eval-meis

Output
------
The output of cluster_identifier is a tab delimited text file with clipped cluster consensus sequences.
The columns are as follows:

|      |                                          |
| ---: | ---------------------------------------- |
| 1.   | Coordinate                               |
| 2.   | Side of read where soft-clipped occurred |
| 3.   | Number of reads in cluster               |
| 4.   | Clipped read consensus                   |
| 5.   | Anchored read consensus                  |
	
The output of SCRAMble-MEIs.R is a tab delimited text file with MEI calls. If no MEIs are present an output file will still be produced with only the header.
The columns are as follows:

|      |                               |                                                                                              |
| ---: | ----------------------------- | -------------------------------------------------------------------------------------------- |
| 1.   | Insertion                     | Coordinate where MEI insertion occurs (zero-based)                                           |
| 2.   | MEI_Family                    | The consensus sequence to which the clipped sequence aligned best                            |
| 3.   | Insertion_Direction           | Whether MEI is on fwd or rev strand relative to bam reference                                |
| 4.   | Clipped_Reads_In_Cluster      | Number of supporting reads in cluster                                                        |
| 5.   | Alignment_Score               | Pairwise alignment score of clipped read consensus to MEI reference sequence                 |
| 6.   | Alignment_Percent_Length      | Percent of clipped read consensus sequence involved in alignment to MEI reference sequence   |
| 7.   | Alignment_Percent_Identity    | Percent identify of alignment of clipped read consensus sequence with MEI reference sequence |
| 8.   | Clipped_Sequence              | Clipped cluster consensus sequences                                                          |
| 9.   | Clipped_Side                  | Left or right, side of read where soft-clipping ocurred                                      |
| 10.  | Start_In_MEI                  | Left-most position of alignment to MEI reference sequence                                    |
| 11.  | Stop_In_MEI                   | Right-most position of alignment to MEI reference sequence                                   |
| 12.  | polyA_Position                | Position of  polyA clipped read cluster if found                                             |
| 13.  | polyA_Seq                     | Clipped cluster consensus sequences of polyA clipped read cluster if found                   |
| 14.  | polyA_SupportingReads         | Number of supporting reads in polyA clipped read cluster if found                            |
| 15.  | TSD                           | Target site duplication sequence if polyA clipped read cluster found                         |
| 16.  | TSD_length                    | Length of target site duplication if polyA clipped read cluster found                        |

R Dependencies
--------------
SCRAMBLE-MEIs.R was developed on R version 3.1.1 (2014-07-10) and uses the following libraries:

    [1] Biostrings_2.34.1   XVector_0.6.0       IRanges_2.0.1
    [4] S4Vectors_0.4.0     BiocGenerics_0.12.1 stringr_0.6.2
    [7] optparse_1.4.4


Disclaimers
-----------
In theory, SCRAMble should work well on any MEI reference fasta sequences, however, it has only been tested on the
sequences provided in `/path/to/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa`.
