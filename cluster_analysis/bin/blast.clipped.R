blast.clipped = function(df, indelScore=80, pctAlign=90, blastRef){
  library(rBLAST)
  
  # align clipped sequences to reference
  bl = blast(db=blastRef)
  
  # create named DNAStringSet object
  seq = df$clipped.consensus
  names(seq) = df$rname_clippedPos_Orientation_ReadSide
  seq = DNAStringSet(seq, use.names=T)
  
  # blast
  results = predict(bl, seq, BLAST_args = "-dust no")
  names(results) = c("query", "subject", "pct_identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  results$qlen = results$qend - results$qstart + 1
  
  ### filter blast results by length and bitscore
  results = results[(results$qlen/results$length)>(pctAlign/100) & results$bitscore > indelScore,]
  
  ##################################
  # compare to mapping of anchored portion of read
  reference.clusters = data.frame(
    chr=df$RNAME,
    reference.start=df$clipped_pos,
    reference.stop=df$clipped_pos,
    rname_clippedPos_Orientation_ReadSide=df$rname_clippedPos_Orientation_ReadSide,
    counts = df$counts,
    strand="+",
    stringsAsFactors = F)
  reference.clusters = merge(reference.clusters, results, by.x="rname_clippedPos_Orientation_ReadSide", by.y="query", all.x=T)

  return(reference.clusters)
}
##################################
##################################
best.hits = function(reference.clusters){
  hits = reference.clusters[!is.na(reference.clusters$subject) & reference.clusters$chr == reference.clusters$subject,]
  hits$hit.len = hits$send-hits$sstart
  
    if(nrow(hits)>0){
      hits$pct_aligned = 100*(hits$qlen )/hits$length
      hits$aligned.strand = sign(hits$send-hits$sstart)
      hits$dist.to.alignment.start = hits$sstart - hits$reference.start
      hits$dist.to.alignment.end   = hits$send - hits$reference.start
      hits$sign.same = sign(hits$dist.to.alignment.end) == sign(hits$dist.to.alignment.start)
      hits = hits[hits$sign.same,] # remove alignments that overlapped anchored portion of read
      hits$start.dist = abs(hits$sstart - hits$reference.start)
      
      ### keep hits with highest score, closest to clipped read
      hits$min.dist = sapply(1:nrow(hits), function(i) min(c(abs(hits$dist.to.alignment.end[i]), abs(hits$dist.to.alignment.start[i]))) )
      hits = hits[order(hits$rname_clippedPos_Orientation_ReadSide, hits$evalue, hits$min.dist),]
      hits = hits[!duplicated(hits$rname_clippedPos_Orientation_ReadSide),]
    }
  return(hits)
}
###########################
