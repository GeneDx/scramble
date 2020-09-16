## USEFUL FUNCTIONS

#########################################
#' @title Reverse Complement
#' @description
#' Given a string or vector of strings, returns the reverse compelement
#'
#' @export
#' @name rc
#' @author Rebecca Torene
#' @param string DNA sequence (A, C, T, G supported)
#' @return reverse complement of string or vector of strings
#' @include IRanges
#' @examples
#' rc('ACGGGACTACTTAACG')
rc = function(string){
  # fxn to take reverse complement of a single string or vector of strings
  library(IRanges)
  r.string = tolower(reverse(string))
  rc.string = gsub("a", "T", r.string)
  rc.string = gsub("t", "A", rc.string)
  rc.string = gsub("c", "G", rc.string)
  rc.string = gsub("g", "C", rc.string)
  return(tolower(rc.string))
}
############################################################
#' @title Complexity Score
#' @description
#' Given a string or vector of strings, returns the chi-square p-value for sequence complexity
#'
#' @export
#' @name complexity.score
#' @author Rebecca Torene
#' @param string DNA sequence (A, C, T, G supported)
#' @return chi-square p-value for sequence complexity with null expectation of random sequence
#' @include Biostrings
#' @examples
#' complexity.score(c('ACGGGACTACTTAACG', 'ACTAAAAAAAAAAAAAAAACCAAAAA'))
complexity.score = function(seq){
  library(Biostrings)
  seq = DNAStringSet(seq) 
  mono = oligonucleotideFrequency(seq, width=1, as.prob=F)
  di   = dinucleotideFrequency(seq, as.prob=F)
  freqs = cbind(mono, di)
  all.mono = rowSums(mono>0)==1
  
  make.pair = function(x){return(c(rep(width(x)/4, 4), rep(width(x)/16, 16)))}
  p = sapply(1:length(seq), function(i) suppressWarnings(chisq.test(cbind(freqs[i,], make.pair(seq[i]))))$p.value)
  p[all.mono] = 0 
  return(p)
}
############################################################
cat("Useful Functions Loaded\n")