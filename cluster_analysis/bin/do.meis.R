## MOBILOME
#Takes the clipped read consensus (forward sequence and its reverse complement) and aligns to the mobilome.
#Alignment quality statistics are returned. Filtering thresholds are used to determine likely MEIs.
library(stringr)
library(Biostrings)
#############################################################
## performs pairwise alignments between vector of sequences and individual reference sequences
## returns data.frame with key stats on alignment
aligner = function(my.fa, ref=mobilome, rc=F){
  results = {}
  for(i in 1:length(ref)){
    alignment=pairwiseAlignment(pattern = tolower(my.fa), subject = ref[[i]], type="local") 
    alignment.df = data.frame(alignment_len=nchar(alignment),
                              alignment_score=score(alignment),
                              percent_identity = pid(alignment), stringsAsFactors = F)
    alignment.df$percent_clipped_read_aligned = (alignment.df$alignment_len * 100) / nchar(my.fa)
    alignment.df$starts = alignment@subject@range@start
    alignment.df$stops = alignment.df$starts + alignment@subject@range@width - 1    

   alignment.df$MEI_Family = names(ref)[i]

    if(i==1){
      results = alignment.df
    }else{
      results = rbind.data.frame(results, alignment.df)
    }
    cat("Done analyzing", names(ref)[i], "\n")
  }
  results$rc = rc
  return(results)
} ## END OF FUNCTION
###################################################################

pick.winners = function(df, pct_read_aligned=75, meiScore=50){
  cols = c("coord","MEI_Family", 
           "alignment_score", "percent_clipped_read_aligned", "percent_identity", "starts", "stops", "rc")
  
  my.winners = df[df$alignment_score >= meiScore & df$percent_clipped_read_aligned >= pct_read_aligned,]
  if(nrow(my.winners)==0){
    winners = data.frame(Insertion=character(), 
                         MEI_Family=character(),
                         Insertion_Direction=character(),
                         Clipped_Reads_In_Cluster=character(),
                         Alignment_Score=character(),
                         Alignment_Percent_Length=character(),
                         Alignment_Percent_Identity=character(),
                         Clipped_Sequence=character(),
                         Clipped_Side=character(),
                         Start_In_MEI=character(),
                         Stop_In_MEI=character(),
                         RNAME= character(),
                         clipped_pos = character(),
                         stringsAsFactors = F)
  }else{
    winners = data.frame(Insertion=my.winners$coord, 
                         MEI_Family=my.winners$MEI_Family,
                         Insertion_Direction=ifelse(!my.winners$rc, "Plus", "Minus"),
                         Clipped_Reads_In_Cluster=my.winners$counts,
                         Alignment_Score=my.winners$alignment_score,
                         Alignment_Percent_Length=my.winners$percent_clipped_read_aligned,
                         Alignment_Percent_Identity=my.winners$percent_identity,
                         Clipped_Sequence=toupper(my.winners$clipped.consensus),
                         Clipped_Side=my.winners$clipped,
                         Start_In_MEI=my.winners$starts,
                         Stop_In_MEI=my.winners$stops,
                         RNAME= my.winners$RNAME,
                         clipped_pos = my.winners$clipped_pos,
                         stringsAsFactors = F)
  }
  return(winners)
} ## END OF FUNCTION
##################################################################
do.meis = function(all,  refs, polyAFrac=0.5, meiScore=50,
                   pctAlign=70, polyAdist=200){
  mobilome = make.mobilome(refs)
  df.all = all[all$counts >= nCluster , ] ## INCLUDE SIMPLE SEQUENCE FOR MEI SEARCH
  df.all$clipped.consensus.rc = rc(df.all$clipped.consensus)
  alignments_fwd = aligner(my.fa=df.all$clipped.consensus, ref=mobilome)
  alignments_rev = aligner(my.fa=df.all$clipped.consensus.rc, ref=mobilome, rc=T)
  df.aligned = rbind.data.frame(data.frame(df.all, alignments_fwd,stringsAsFactors = F, row.names=NULL),
                                data.frame(df.all, alignments_rev, stringsAsFactors = F, row.names=NULL) )

  ## Make pretty output table
  df.aligned$coord = paste(df.all$RNAME, df.all$clipped_pos, sep=":")
  winners = as.data.frame(pick.winners(df.aligned,  pct_read_aligned=pctAlign, meiScore=meiScore) )
  cat("Sample had", nrow(winners), "MEI(s)\n")
  
  
  ###################################################
  ## IDENTIFY POLY-A TAIL AND TSD
  
  if(nrow(winners)==0){
    winners = data.frame(winners,
                         polyA_Position=character(), polyA_Seq=character(),
                         polyA_SupportingReads=character(), TSD=character(),
                         TSD_length=character())
  }else{
    winners$polyA_Position = "None Found"
    winners$polyA_Seq = "None Found"
    winners$polyA_SupportingReads = "None Found"
    winners$TSD = "None Found"
    winners$TSD_length = "NA"
  }
  if(nrow(winners)>0){
    for(i in 1:nrow(winners)){
      if(is.na(winners$Insertion_Direction[i]) | winners$Insertion_Direction[i] == "Plus"){
        # then look for left clipped with poly-A upstream
        polyA = all[all$clipped=="left" & all$RNAME == winners$RNAME[i] & all$clipped_pos >= winners$clipped_pos[i]-polyAdist & all$clipped_pos <= winners$clipped_pos[i],]
        if(nrow(polyA)==0){ next }
        polyA$A_count = str_count(polyA$clipped.consensus, "a")
        polyA$A_frac = polyA$A_count/nchar(polyA$clipped.consensus)
        if( any(!is.na(polyA$A_frac) & polyA$A_frac >= polyAFrac)){
          # grab the one closest to the MEI
          polyA=polyA[order(polyA$clipped_pos, decreasing=T),]
          polyA = polyA[polyA$A_frac>=polyAFrac,][1,]
          winners$TSD[i] = toupper(substr(polyA$anchored.consensus, start=1, stop=(winners$clipped_pos[i] - polyA$clipped_pos ) ))
          winners$TSD_length[i] = winners$clipped_pos[i] - polyA$clipped_pos
          if((winners$clipped_pos[i] - polyA$clipped_pos ) > nchar(polyA$anchored.consensus)){
            winners$TSD[i] = paste0(winners$TSD[i], "<TRUNCATED>",  collapse="")
          }
        }else{next()}
      }else{
        # then look for right clipped with poly-T downstream
        polyA = all[all$clipped=="right" & all$RNAME == winners$RNAME[i] & all$clipped_pos <= winners$clipped_pos[i]+polyAdist & all$clipped_pos >= winners$clipped_pos[i],]
        if(nrow(polyA)==0){ next }
        polyA$A_count = str_count(polyA$clipped.consensus, "t")
        polyA$A_frac = polyA$A_count/nchar(polyA$clipped.consensus)
        if(any(!is.na(polyA$A_frac) & polyA$A_frac >= polyAFrac)){
          # grab the one closest to the MEI
          polyA=polyA[order(polyA$clipped_pos, decreasing=F),]
          polyA = polyA[polyA$A_frac>=polyAFrac,][1,]
          winners$TSD_length[i] = polyA$clipped_pos - winners$clipped_pos[i]
          if(1 > (nchar(polyA$anchored.consensus) - (polyA$clipped_pos - winners$clipped_pos[i]))){
            winners$TSD[i] = paste0("<TRUNCATED>",toupper(polyA$anchored.consensus), collapse="")
          }else{
            winners$TSD[i] = toupper(substr(polyA$anchored.consensus, start= (nchar(polyA$anchored.consensus) - (polyA$clipped_pos - winners$clipped_pos[i]) + 1), stop=nchar(polyA$anchored.consensus)) )
          }
        } else{next()}
      }
      winners$polyA_Position[i] = polyA$clipped_pos
      winners$polyA_Seq[i] = toupper(polyA$clipped.consensus)
      winners$polyA_SupportingReads[i] = polyA$counts
      
    }
  }
  
  
  winners = winners[,-(match(c("RNAME", "clipped_pos"), names(winners)) )]
  winners = winners[order(winners$Insertion, winners$Alignment_Score, decreasing=T),]
  winners = winners[!duplicated(winners$Insertion),]
  return(winners)
}  
##################################################################
## MAKE MOBILOME
# convert fasta reference MEI sequences into mobilome R object
make.mobilome = function(mei.refs){
  refs = read.delim(mei.refs, as.is=T, header=F)
  mobilome = list()
  indices = which(grepl("^>", refs[,1]))
  for(i in 1:length(indices)){
    seq = paste0(refs[(indices[i]+1):ifelse(i<length(indices),(indices[i+1]-1),nrow(refs)),1], collapse="")
    mobilome = c(mobilome, list(tolower(seq)))
    names(mobilome)[i] = gsub(">", "", as.character(refs[indices[i],1]))
  }
  return(mobilome)
}
