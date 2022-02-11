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
get_refs = function(fa, chrom, start, end){
  if (missing(fa) | missing(chrom) | missing(start) | missing(end)) return('N')
  if (! chrom %in% names(fa)) return('N')
  fa = fa[chrom]
  seq = subseq(fa, start=start, end=end)
  return(as.vector(seq))
}
##############################
make.vcf.header = function(fa, blastRef=None){
  if (missing(fa)) return(NULL)
  contigs = names(fa)
  header = c('##fileformat=VCFv4.2',
             paste('##reference=', blastRef, sep=''),
             paste('##contig=<ID=', contigs, '>', sep=""),
             '##FILTER=<ID=PASS,Description="All filters passed">',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
             '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
             '##INFO=<ID=END,Number=.,Type=Integer,Description="End position for structural variants">',
             '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
             '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
             '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
             '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">'
             )
  return(header)
}

##############################
write.scramble.vcf = function(winners, fa, meis=F){

    # return empty fixed data when no variants found
    if(nrow(winners) == 0){
        fixed = data.frame('#CHROM' = character(),
               POS = character(),
               ID = character(),
               REF = character(),
               ALT = character(),
               QUAL = character(),
               FILTER = character(),
               INFO = character(),
               check.names = F)
        return(fixed)
    }

  #argument checks
  if (is.null(winners)) return(NULL)

  if(!meis){
    fixed = data.frame('#CHROM' = winners$CONTIG,
                       POS = winners$DEL.START,
                       ID = 'DEL',
                       QUAL = sapply(1:nrow(winners), function(i) get_score(winners$SCORE.RIGHT.ALIGNMENT[i], winners$SCORE.LEFT.ALIGNMENT[i])),
                       FILTER = 'PASS',
                       REF = sapply(1:nrow(winners), function(i) get_refs(fa, winners$CONTIG[i], winners$DEL.START[i], winners$DEL.END[i] + 1)),
                       svtype = 'DEL',
                       stringsAsFactors = F, check.names = F)
    
    fixed$ALT = str_sub(fixed$REF, start=-1)
    fixed$svlen = nchar(fixed$REF)
    fixed$end = fixed$POS + fixed$svlen
    fixed$INFO = paste('SVTYPE=', fixed$svtype, ';', 'SVLEN=', fixed$svlen, ';', 'END=', fixed$end, sep='')
  } else {
    fixed = data.frame('#CHROM' =  gsub("(.*):(\\d*)$", "\\1", winners$Insertion),
                       POS = as.integer(gsub("(.*):(\\d*)$", "\\2", winners$Insertion)),
                       ID = 'INS:ME',
                       FILTER = 'PASS',
                       ALT = paste('<INS:ME:', toupper(winners$MEI_Family), '>', sep=''),
                       QUAL = winners$Alignment_Score,
                       name = paste(winners$Insertion, toupper(winners$MEI_Family), winners$Insertion_Direction, sep="_"),
                       polarity = ifelse(winners$Insertion_Direction == 'Plus', "+", "-"),
                       stringsAsFactors = F, check.names = F)
    fixed$start = fixed$POS
    fixed$end = fixed$POS+1
    fixed$INFO = paste('MEINFO', paste(fixed$name, fixed$start, fixed$end, fixed$polarity, sep=','), sep='=')
    fixed$REF = sapply(1:nrow(fixed), function(i) get_refs(fa, fixed[i, '#CHROM'], fixed$POS[i], fixed$POS[i]))
  }   

  vcf.cols = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
  return(fixed[,vcf.cols])

}
