#!/usr/bin/env Rscript
##############################################
# This script performs SCRAMble SV analysis  #
# Developed by R.Torene                      #
                                             #
# parses command line options                #
# reads in cluster files                     #
# performs pre-processing of data            #
# analyzes clusters for MEIs and Deletions   #
##############################################
library(optparse)
parser <- OptionParser()
## analysis options
parser <- add_option(parser, c("-o", "--out-name"), type="character",
                     dest="out.file.prefix", help="full path to output file including prefix")
parser <- add_option(parser, c("-c", "--cluster-file"), type="character",
                     dest="cluster.file",  help="full path to cluster file")
parser <- add_option(parser, c("-n", "--nCluster"), type="integer",
                     dest="nCluster", default=5, help="min cluster size to analyze [default %default]")
parser <- add_option(parser, c("--indel-score"),  type="integer", default=80,
                     dest="indel.score", help="min indel alignment score to call [default %default]")
parser <- add_option(parser, c("-m", "--mei-score"),  type="integer", default=50,
                     dest="mei.score", help="min mei alignment score to call [default %default]")
parser <- add_option(parser, c("--min-del-len"),  type="integer", default=50,
                     dest="min.del.len", help="minimum deletion length to call [default %default]")
parser <- add_option(parser, c("--poly-a-frac"),  type="numeric", default=0.75,
                     dest="polyAFrac", help="fraction of clipped length for calling polyA tail in MEIs [default %default]")
parser <- add_option(parser, c("--pct-align"),  type="integer", default=90,
                     dest="pctAlign", help="percent alignment of clipped read for calling deletionss [default %default]")
parser <- add_option(parser, c("--poly-a-dist"),  type="integer", default=100,
                     dest="polyAdist", help="how far from MEI to look for polyA tail [default %default]")
parser <- add_option(parser, c("-i", "--install-dir"), type="character",
                     default="~/scramble/cluster_analysis/bin/",
                     dest="install.dir", help="resource directory [default %default]")
parser <- add_option(parser, c("--mei-refs"), type="character",
                     default="~/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa",
                     dest="mei.ref.file", help="full path to MEI reference file (fasta format) [default %default]")
parser <- add_option(parser, c("-r", "--ref"), type="character",
                     default=NULL, dest="ref", help="reference file (fasta format) [default %default]")
parser <- add_option(parser, c("--no-vcf"), action="store_true", default=FALSE,
                     type="logical", help="evaluate meis")
## what to evaluate
parser <- add_option(parser, c("--eval-meis"), action="store_true", default=FALSE,
                     type="logical", help="evaluate meis")
parser <- add_option(parser, c("--eval-dels"), action="store_true", default=FALSE,
                     type="logical", help="evaluate deletions")
opt = parse_args(parser)
###############################
## init globals from arguments
outFilePrefix   = opt[["out.file.prefix"]]                 
clusterFile     = opt[["cluster.file"]]                 
nCluster        = opt[["nCluster"]]  
indelScore      = opt[['indel.score']]
meiScore        = opt[["mei.score"]]  
minDelLen       = opt[['min.del.len']]
mei.refs        = opt[["mei.ref.file"]]         
INSTALL.DIR     = opt[["install.dir"]]   
polyAFrac       = opt[["polyAFrac"]]  
pctAlign        = opt[["pctAlign"]]  
polyAdist       = opt[["polyAdist"]] 
blastRef        = opt[["ref"]]
no.vcf          = opt[["no-vcf"]]

# what to evaluate
meis        = opt[["eval-meis"]]
deletions   = opt[["eval-dels"]]
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
###############################
## if no eval options provided log and exit
if(!meis & !deletions){
    cat('No structural variants to evaluate. Please use flags --eval-meis and --eval-dels to indicate analysis.')
    stop()
}
###############################
## READ IN DATA AND PRE-PROCESS
all = read.delim(clusterFile, as.is=T, header=F,
                 col.names=c("coord", "clipped", "counts", "clipped.consensus", "anchored.consensus"))
all$clipped_pos = as.integer(gsub(".*:", "", all$coord) )
all$RNAME = gsub("(.*):\\d*$", "\\1", all$coord)
all$rname_clippedPos_Orientation_ReadSide = paste(all$RNAME, all$clipped_pos, all$clipped, sep="_")
##############################
## MEIS
if(meis){
  source("do.meis.R")
  mei.winners=do.meis(all=all, meiScore=meiScore, refs=mei.refs,
                  polyAFrac=polyAFrac, pctAlign=pctAlign, polyAdist=polyAdist)
  cat("Done analyzing MEIs\n")
  write.table(mei.winners, paste(outFilePrefix, "_MEIs.txt", sep=""), row.names=F, quote=F, sep="\t")  
}
##############################
## DELETIONS
if(deletions){
  source("do.dels.R")
  del.winners = do.dels(all, indelScore, minDelLen, pctAlign, blastRef)
  write.table(del.winners, paste(outFilePrefix, "_PredictedDeletions.txt", sep=""), sep="\t", quote=F, row.names=F)
  cat("Done analyzing deletions\n")
}
##############################
## WRITE VCF
if(!no.vcf){
  if(is.null(blastRef)){
    warning("A reference .fa file is required for writing results to VCF")
  }else{
    source('make.vcf.R')
    message("Writing VCF file to ", paste0(outFilePrefix, ".vcf"), "...")

    #Load reference fasta
    suppressMessages(require(Rsamtools))
    fa = getSeq(open(FaFile(blastRef)))
    #Contig ID = FASTA ID until first space character
    names(fa) = gsub(" .*", "", names(fa))

    write.table(make.vcf.header(fa, blastRef), paste0(outFilePrefix, ".vcf"), row.names=F, col.names=F, quote=F)

    # get mei results to fixed data frame
    mei.fixed = NULL
    if (meis) mei.fixed = write.scramble.vcf(mei.winners, fa, T)
    
    # get del results to fixed data frame
    del.fixed = NULL
    if (deletions) del.fixed = write.scramble.vcf(del.winners, fa, F)
    
    fixed = rbind.data.frame(mei.fixed, del.fixed)
    fixed = fixed[order(fixed[,'#CHROM'], fixed$POS), ]
    suppressWarnings(write.table(fixed, paste0(outFilePrefix, ".vcf"), row.names=F, col.names=T, quote=F, append=T, sep='\t'))
    message("Done.")
  }
}




