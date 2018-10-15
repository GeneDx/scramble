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
#########################################
cat("Useful Functions Loaded\n")