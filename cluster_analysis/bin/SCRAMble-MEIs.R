#!/usr/bin/env Rscript
##############################################
# This script performs SCRAMble MEI analysis #
# Developed by R.Torene                      #
                                             #
# parses command line options                #
# reads in cluster files                     #
# performs pre-processing of data            #
# analyzes clusters for MEIs                 #
##############################################
library(optparse)
parser <- OptionParser()
## analysis options
parser <- add_option(parser, c("-o", "--out-name"), type="character",
                     dest="out.file", help="full path to output file")
parser <- add_option(parser, c("-c", "--cluster-file"), type="character",
                     dest="cluster.file",  help="full path to cluster file")
parser <- add_option(parser, c("-n", "--nCluster"), type="integer",
                     dest="nCluster", default=5, help="min cluster size to analyze [default %default]")
parser <- add_option(parser, c("-m", "--min-score"),  type="integer", default=50,
                     dest="mei.score", help="min MEI alignment score to make a call [default %default]")
parser <- add_option(parser, c("--pct-align"),  type="integer", default=70,
                     dest="pctAlign", help="min percent of clipped sequence length aligned to MEI [default %default]")
parser <- add_option(parser, c("--poly-a-frac"),  type="numeric", default=0.5,
                     dest="polyAFrac", help="min fraction of clipped seqeunce length made of As to call polyA tail [default %default]")
parser <- add_option(parser, c("--poly-a-dist"),  type="integer", default=200,
                     dest="polyAdist", help="max distance from MEI to search for polyA tail [default %default]")
parser <- add_option(parser, c("-i", "--install-dir"), type="character",
                     default="~/scramble/cluster_analysis/bin/",
                     dest="install.dir", help="resource directory [default %default]")
parser <- add_option(parser, c("--mei-refs"), type="character",
                     default="~/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa",
                     dest="mei.ref.file", help="full path to MEI reference file (fasta format) [default %default]")
opt = parse_args(parser)
###############################
## init globals from arguments
outFile         = opt[["out.file"]]                 
clusterFile     = opt[["cluster.file"]]                 
nCluster        = opt[["nCluster"]]                 
meiScore        = opt[["mei.score"]]     
mei.refs        = opt[["mei.ref.file"]]         
INSTALL.DIR     = opt[["install.dir"]]   
polyAFrac       = opt[["polyAFrac"]]  
pctAlign        = opt[["pctAlign"]]  
polyAdist       = opt[["polyAdist"]]  
###############################
cat("Running sample:", clusterFile, "\n")
setwd(INSTALL.DIR)
objects = ls()
objects = objects[-match(c("opt", "parser"),objects)]
# print(objects)
cat("Running scramble with options:\n")
for(i in 1:length(objects)){
  cat(objects[i],":", get(objects[i]), "\n")
}
source("usefulFunctions.R")
source("do.meis.R")
# library(tools)
###############################
## READ IN DATA AND PRE-PROCESS
all = read.delim(clusterFile, as.is=T, header=F,
                 col.names=c("coord", "clipped", "counts", "clipped.consensus", "anchored.consensus"))
all$clipped_pos = as.integer(gsub(".*:", "", all$coord) )
all$RNAME = gsub(":.*", "", all$coord)
all$rname_clippedPos_Orientation_ReadSide = paste(all$RNAME, all$clipped_pos, all$clipped, sep="_")
##############################
## MEIS
winners=do.meis(all=all, meiScore=meiScore, refs=mei.refs,
                polyAFrac=polyAFrac, pctAlign=pctAlign, polyAdist=polyAdist)
cat("Done analyzing MEIs\n")
write.table(winners, outFile, row.names=F, quote=F, sep="\t")  
# print(sessionInfo())
