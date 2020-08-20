##############################
# convert to VCF
##############################
get_score = function(right_score, left_score){
  if(is.na(right_score)){
    return(left_score)
  }else if(is.na(left_score)){
    return(right_score)
  }else{
    return(mean(c(left_score, right_score)))
  }
}
##############################
get_refs = function(blastRef, start, end){
  fa = getSeq(open(FaFile(blastRef)))
  seq = subseq(fa, start=start, end=end)
  return(as.vector(seq))
}
##############################
write.scramble.vcf = function(outFile, del.winners){
  ## TODO: write VCF for SCRAMble MEI calls too
  fixed = data.frame('#CHROM' = del.winners$CONTIG,
                     POS = del.winners$DEL.START,
                     ID = 'DEL',
                     REF = sapply(1:nrow(del.winners), function(i) get_refs(blastRef, del.winners$DEL.START[i], del.winners$DEL.END[i] + 1)),
                     QUAL = sapply(1:nrow(del.winners), function(i) get_score(del.winners$SCORE.RIGHT.ALIGNMENT[i], del.winners$SCORE.LEFT.ALIGNMENT[i])),
                     FILTER = 'PASS',
                     reads = sapply(1:nrow(del.winners), function(i) get_score(del.winners$RIGHT.CLUSTER.COUNTS[i], del.winners$LEFT.CLUSTER.COUNTS[i])),
                     svtype = 'DEL',
                     stringsAsFactors = F, check.names = F)
  fixed$ALT = str_sub(fixed$REF, start=-1)
  fixed$svlen = nchar(fixed$REF)
  fixed$end = fixed$POS + fixed$svlen
  fixed$INFO = paste('SVTYPE=', fixed$svtype, ';', 'SVLEN=', fixed$svlen, ';', 'END=', fixed$end, sep='')
  contigs = paste('contig=<ID=', unique(fixed[,'#CHROM']), '>', sep='')
  fixed = fixed[order(fixed[,'#CHROM'], fixed$POS), ]
  
  header = c('##fileformat=VCFv4.2',
             paste('##reference=', blastRef, sep=''),
             contigs,
             '##FILTER=<ID=PASS,Description="All filters passed">',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
             '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
             '##INFO=<ID=END,Number=.,Type=Integer,Description="End position for structural variants">'
             )
  vcf.cols = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
  write.table(header, paste(outFilePrefix, ".vcf", sep=""), row.names=F, col.names=F, quote=F)
  suppressWarnings(write.table(fixed[,vcf.cols], paste(outFilePrefix, ".vcf", sep=""), row.names=F, col.names=T, quote=F, append=T, sep='\t'))
}