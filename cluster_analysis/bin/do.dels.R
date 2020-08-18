## FIND DELETIONS
del.finder = function(dat){
  print(paste("Two-End-Deletions: Working on contig", unique(dat$RNAME)))
  # find right clipped clusters upstream of left clipped clusters
  dat = dat[order(dat$clipped_pos),]
  dat=dat[!duplicated(dat$rname_clippedPos_Orientation_ReadSide),]
  rights = which(dat$clipped=="right")
  lefts = which(dat$clipped=="left")
  dat$clipped.len = nchar(dat$clipped.consensus)
  dat$anchored.len = nchar(dat$anchored.consensus)
  my.rights=rights[(rights+1) %in% lefts]
  
  # initiate alignments object:
  alignments=data.frame(CONTIG = character(),
                        DEL.START = character(),
                        DEL.END = character(),
                        REF.ANCHOR.BASE = character(),
                        DEL.LENGTH = character(),
                        RIGHT.CLUSTER=character(),
                        RIGHT.CLUSTER.COUNTS = character(),
                        LEFT.CLUSTER=character(),
                        LEFT.CLUSTER.COUNTS = character(),
                        LEN.RIGHT.ALIGNMENT=character(),
                        SCORE.RIGHT.ALIGNMENT=character(),
                        PCT.COV.RIGHT.ALIGNMENT=character(),
                        PCT.IDENTITY.RIGHT.ALIGNMENT=character(),
                        LEN.LEFT.ALIGNMENT=character(),
                        SCORE.LEFT.ALIGNMENT=character(),
                        PCT.COV.LEFT.ALIGNMENT=character(),
                        PCT.IDENTITY.LEFT.ALIGNMENT=character(),
                        INS.SIZE = character(),
                        INS.SEQ = character(),
                        RIGHT.CLIPPED.SEQ = character(),
                        LEFT.CLIPPED.SEQ = character())
  if(length(my.rights)==0){return(alignments)}
  for (i in my.rights){
    right.clipped = dat$clipped.consensus[i]
    down.cluster = paste0(dat$clipped.consensus[i + 1], dat$anchored.consensus[i + 1], collapse="")
    my.right.alignment=pairwiseAlignment(pattern = right.clipped, subject = down.cluster, type="local")
    left.clipped = dat$clipped.consensus[i + 1]
    up.cluster = paste0(dat$anchored.consensus[i], dat$clipped.consensus[i], collapse="")
    my.left.alignment=pairwiseAlignment(pattern = left.clipped, subject = up.cluster, type="local")
    
    #### DETERMINE INSERTION SIZE
    ## finds deletions with short insertions
    cov=as.vector(coverage(my.right.alignment))
    # if first 1 is at position > length downstream clipped, then no inserted seq
    if(min(which(cov==1))>dat$clipped.len[i+1]){
      INS.SIZE=0
      INS.SEQ=NA
    }else if(max(which(cov==1))<=dat$clipped.len[i+1]){
      # if last 1 is at position < length downstream clipped, then long insertion
      # ins is length of alignment + unaligned portions of clipped sequences
      INS.SIZE =  dat$clipped.len[i] + (length(cov) - dat$anchored.len[i+1] - nchar(my.right.alignment))
      INS.SEQ = toupper(paste0(dat$clipped.consensus[i], substr(dat$clipped.consensus[i+1], start =  nchar(my.right.alignment)+1, stop=dat$clipped.len[i+1]), collapse=""))
    }else{
      # if last 1 is at position > length downstream clipped, then short insertion
      # ins is length of alignment        
      INS.SIZE = sum(cov[1:dat$clipped.len[i+1]])
      INS.SEQ = toupper(substr(dat$clipped.consensus[i+1], start = min(which(cov==1)), stop=dat$clipped.len[i+1]))
    }
    
    
    alignments = rbind.data.frame(alignments, data.frame(
      CONTIG = dat$RNAME[i],
      DEL.START = dat$clipped_pos[i],
      DEL.END = dat$clipped_pos[i + 1],
      REF.ANCHOR.BASE = toupper(substr(dat$anchored.consensus[i], start=dat$anchored.len[i], stop=dat$anchored.len[i])),
      DEL.LENGTH = dat$clipped_pos[i + 1] - dat$clipped_pos[i ] + 1,
      RIGHT.CLUSTER=dat$rname_clippedPos_Orientation_ReadSide[i ],
      RIGHT.CLUSTER.COUNTS = as.integer(dat$counts[i ]),
      LEFT.CLUSTER=dat$rname_clippedPos_Orientation_ReadSide[i+1],
      LEFT.CLUSTER.COUNTS = as.integer(dat$counts[i + 1]),
      LEN.RIGHT.ALIGNMENT=nchar(my.right.alignment),
      SCORE.RIGHT.ALIGNMENT=score(my.right.alignment),
      PCT.COV.RIGHT.ALIGNMENT=(nchar(my.right.alignment) * 100) / nchar(right.clipped),
      PCT.IDENTITY.RIGHT.ALIGNMENT=pid(my.right.alignment),
      LEN.LEFT.ALIGNMENT=nchar(my.left.alignment),
      SCORE.LEFT.ALIGNMENT=score(my.left.alignment),
      PCT.COV.LEFT.ALIGNMENT=(nchar(my.left.alignment) * 100) / nchar(left.clipped),
      PCT.IDENTITY.LEFT.ALIGNMENT=pid(my.left.alignment),
      INS.SIZE = INS.SIZE,
      INS.SEQ = INS.SEQ,
      RIGHT.CLIPPED.SEQ = toupper(right.clipped),
      LEFT.CLIPPED.SEQ = toupper(left.clipped), stringsAsFactors = F  ))
  }
  return(alignments)
}

two.end.dels = function(df, minDelLen = 50, indelScore){
  #Look for right-clipped reads that are upstream of left-clipped reads in the same contig. Align the upstream, right-clipped sequence to the downstream, left-clipped cluster sequence. Conversely, align the downstream, left-clipped sequence to the upstream, right-clipped cluster sequence. High scoring alignements mean a simple deletion. 
  #In looking for deletions with short insertions, the alignment will initiate within the clipped portion of the other cluster.
  if(nrow(df)==0){
    byContig = list(df)
  }else{
    byContig =  split(df, df$RNAME)
  }
    
  x = lapply(byContig, del.finder)
  z = do.call(rbind, lapply(x, data.frame, stringsAsFactors=FALSE))
  cols=c("DEL.LENGTH",  "RIGHT.CLUSTER.COUNTS", "LEFT.CLUSTER.COUNTS", "LEN.RIGHT.ALIGNMENT", "SCORE.RIGHT.ALIGNMENT", "PCT.COV.RIGHT.ALIGNMENT", "PCT.IDENTITY.RIGHT.ALIGNMENT", "LEN.LEFT.ALIGNMENT", "SCORE.LEFT.ALIGNMENT", "PCT.COV.LEFT.ALIGNMENT", "PCT.IDENTITY.LEFT.ALIGNMENT", "INS.SIZE")
  
  if(nrow(z)==0){
    delWinners=z
  }else{
    z[,cols] = apply(z[,cols], 2, function(x) as.numeric(as.character(x)) )
    delWinners = z[!is.na(z$CONTIG) & z$SCORE.LEFT.ALIGNMENT > indelScore & z$SCORE.RIGHT.ALIGNMENT > indelScore & z$DEL.LENGTH >= minDelLen,]
  }
  
  return(delWinners)
}

## ONE END DELETIONS USING BLAST
#BLAST clipped sequence to reference. Keep highest scoring alignments that are closest to the originating read.
#################################
one.end.dels=function(hits, minDelLen=50){

  if(nrow(hits)>0){
    hits$downstream = sapply(1:nrow(hits), function(i) ifelse(hits$dist.to.alignment.start[i]>0, T, F) )
  }else{
    hits$downstream = as.logical(NULL)
  }
  
  
  hits$clipped = gsub("^.*_\\d.*_", "", hits$rname_clippedPos_Orientation_ReadSide)
  hits$pos = as.integer(gsub("(^.*_)(\\d.*)(_.*$)", "\\2", hits$rname_clippedPos_Orientation_ReadSide))
  hits$contig = gsub("_.*", "", hits$rname_clippedPos_Orientation_ReadSide)
  # if R clipped and gene on + strand, then keep downstream alignments on + strand
  # if R clipped and gene on - strand, then keep upstream alignments on - strand
  # if L clipped and gene on + strand, then keep upstream alignments on + strand
  # if L clipped and gene on - strand, then keep downstream alignments on - strand
  hits = hits[(hits$clipped=="right" & hits$strand=="+" & hits$downstream & hits$aligned.strand>0) |
                (hits$clipped=="right" & hits$strand=="-" & !hits$downstream & hits$aligned.strand<0) |
                (hits$clipped=="left" & hits$strand=="+" & !hits$downstream & hits$aligned.strand>0) |
                (hits$clipped=="left" & hits$strand=="-" & hits$downstream & hits$aligned.strand<0),]
  
  # add back info about clipped clusters
  hits$counts = df$counts[match(hits$rname_clippedPos_Orientation_ReadSide, df$rname_clippedPos_Orientation_ReadSide)]
  hits$clipped.consensus = df$clipped.consensus[match(hits$rname_clippedPos_Orientation_ReadSide, df$rname_clippedPos_Orientation_ReadSide)]
  hits$anchored.consensus = df$anchored.consensus[match(hits$rname_clippedPos_Orientation_ReadSide, df$rname_clippedPos_Orientation_ReadSide)]

  ### write results
  alignments = data.frame(CONTIG = character(),
                          DEL.START = character(),
                          DEL.END = character(),
                          REF.ANCHOR.BASE = character(),
                          DEL.LENGTH = character(),
                          RIGHT.CLUSTER = character(),
                          RIGHT.CLUSTER.COUNTS = character(),
                          LEFT.CLUSTER = character(),
                          LEFT.CLUSTER.COUNTS = character(),
                          LEN.RIGHT.ALIGNMENT=character(),
                          SCORE.RIGHT.ALIGNMENT=character(),
                          PCT.COV.RIGHT.ALIGNMENT=character(),
                          PCT.IDENTITY.RIGHT.ALIGNMENT=character(),
                          LEN.LEFT.ALIGNMENT=character(),
                          SCORE.LEFT.ALIGNMENT=character(),
                          PCT.COV.LEFT.ALIGNMENT=character(),
                          PCT.IDENTITY.LEFT.ALIGNMENT=character(),
                          INS.SIZE = character(),
                          INS.SEQ = character(),
                          RIGHT.CLIPPED.SEQ = character(),
                          LEFT.CLIPPED.SEQ = character(), stringsAsFactors=F )
  
  if(nrow(hits)>0){
    for (i in 1:nrow(hits)){
      del.len = min(c(abs(hits$dist.to.alignment.start[i] ), abs(hits$dist.to.alignment.end[i])))
      
      if(hits$downstream[i]){
        del.start = hits$reference.start[i]
        del.stop = del.start + del.len - 1
      }else{
        del.stop = hits$reference.start[i]
        del.start = del.stop - del.len + 1  
      }
      
      alignments = rbind.data.frame(alignments, data.frame(
        CONTIG = hits$chr[i],
        DEL.START = del.start,
        DEL.END = del.stop,
        REF.ANCHOR.BASE = NA,
        DEL.LENGTH = del.len,
        RIGHT.CLUSTER = ifelse(hits$clipped[i]=="right", hits$rname_clippedPos_Orientation_ReadSide[i], NA),
        RIGHT.CLUSTER.COUNTS = ifelse(hits$clipped[i]=="right", hits$counts[i], NA),
        LEFT.CLUSTER = ifelse(hits$clipped[i]=="left", hits$rname_clippedPos_Orientation_ReadSide[i], NA),
        LEFT.CLUSTER.COUNTS = ifelse(hits$clipped[i]=="left", hits$counts[i], NA),
        LEN.RIGHT.ALIGNMENT=ifelse(hits$clipped[i]=="right", hits$qlen[i], NA),
        SCORE.RIGHT.ALIGNMENT=ifelse(hits$clipped[i]=="right", hits$bitscore[i],NA),
        PCT.COV.RIGHT.ALIGNMENT=ifelse(hits$clipped[i]=="right",hits$pct_aligned[i],NA),
        PCT.IDENTITY.RIGHT.ALIGNMENT=ifelse(hits$clipped[i]=="right", hits$pct_identity[i],NA),
        LEN.LEFT.ALIGNMENT=ifelse(hits$clipped[i]=="left", hits$qlen[i], NA),
        SCORE.LEFT.ALIGNMENT=ifelse(hits$clipped[i]=="left", hits$bitscore[i],NA),
        PCT.COV.LEFT.ALIGNMENT=ifelse(hits$clipped[i]=="left",hits$pct_aligned[i], NA),
        PCT.IDENTITY.LEFT.ALIGNMENT=ifelse(hits$clipped[i]=="left", hits$pct_identity[i],NA),
        INS.SIZE = NA,
        INS.SEQ = NA,
        RIGHT.CLIPPED.SEQ = toupper(ifelse(hits$clipped[i]=="right", toupper(hits$clipped.consensus[i]), NA)),
        LEFT.CLIPPED.SEQ = toupper(ifelse(hits$clipped[i]=="left", toupper(hits$clipped.consensus[i]), NA)), stringsAsFactors = F))
    }
  }

  alignments$DEL.LENGTH = as.integer(alignments$DEL.LENGTH)
  alignments = alignments[alignments$DEL.LENGTH >= minDelLen,]
  return(alignments)
}
###########################################################
do.dels = function(df,indelScore, minDelLen=50, pctAlign=90, blastRef){
  # REMOVE SIMPLE CLIPPED READ CLUSTERS AND ALIGN TO reference
  if(nrow(df)==0){ # account for 0 clipped clusters
    df$complexity.p = numeric(0)
  }else{
    df$complexity.p = complexity.score(df$clipped.consensus)
  }
  before = nrow(df)
  df = df[df$complexity.p >= 0.00001,]
  after = nrow(df)
  cat(before-after, "clusters out of", before, "were removed due to simple sequence\n")
  source("blast.clipped.R")
  
  aligned.clusters=blast.clipped(df, indelScore=indelScore, pctAlign=pctAlign, blastRef)

  cat("Number of alignments meeting thresholds:", nrow(aligned.clusters), "\n")
  hits = best.hits(aligned.clusters)
  cat("Number of best alignments:", nrow(hits), "\n")

  two.end.del.winners=two.end.dels(df, minDelLen, indelScore)
  one.end.del.winners=one.end.dels(hits, minDelLen=minDelLen)
  
  # keep only one-end deletions that are NOT included in the two-del deletions
  all.delWinners = rbind.data.frame(two.end.del.winners, one.end.del.winners[!one.end.del.winners$RIGHT.CLUSTER %in% two.end.del.winners$RIGHT.CLUSTER & !one.end.del.winners$LEFT.CLUSTER %in% two.end.del.winners$LEFT.CLUSTER,])
  cat("Sample had", nrow(all.delWinners), "deletions\n")
  return(all.delWinners)
}
###########################################################


