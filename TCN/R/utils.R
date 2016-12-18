#General utility functions.

#' Standardise chromosome names
#'
#' Tuned for human, standardise chromosome names.
#'
#' @export
#' @param chrs The chromosomes to fix up.
#' @return The cleaned up chromosome names.
fix_chr = function(chrs) {
  chrs = as.character(chrs)
  chrs = gsub('^(chr|chrom|chromosome)','',chrs)
  chrs[chrs=='23'] = 'X'
  chrs[chrs=='24'] = 'Y'
  chrs[chrs=='25'] = 'M'
  chrs[chrs=='MT'] = 'M'
  return(chrs)
}

#' Add a prefix to chromosomes
#'
#' Strips any existing prefix and appends a prefix of your choice.
#'
#' @export
#' @param chrs The chromosome names to prefix.
#' @param prefix The prefix to prepend.
#' @return The processed chromosome names.
prefix_chr = function(chrs,prefix='chr') {
  chrs = fix_chr(chrs)
  if(length(chrs)==0){
    return(chrs)
  }else{
    return(paste0(prefix,chrs))
  }
}

#' Make chromosome names into integers
#'
#' Converts chromosome names into integers, useful for sorting.
#'
#' @param chrs The chromosomes to convert.
#' @return The converted chromosomes as a numeric vector.
num_chr = function(chrs) {
  chrs = fix_chr(chrs)
  chrs[chrs=='X'] = '23'
  chrs[chrs=='Y'] = '24'
  chrs[chrs=='M'] = '25'
  return(as.numeric(chrs))
}

#' Replace variables in string
#'
#' An extended version of sprintf that allows named variables to be replaced within a string.
#'
#' @export
#' @param string The string to replace variables in.
#' @param ... What to replace with what.
#' @return The parsed string.
#' @examples
#' varString("It's %day% today in the month of %month%",day="Monday",month="Smarch")
varString = function(string,...) {
  out = string
  vars = list(...)
  for(i in seq_along(vars)){
    out = gsub(paste0('%',names(vars)[i],'%'),vars[[i]],out)
  }
  return(out)
}

#' An enhancement to the standard table function 
#'
#' Pads out standard calls to table with extra levels not present in the data without the need to force input data to be factors.
#'
#' @export
#' @param ... The thing to be counted, passed to \code{table}.
#' @param levels These are entries that are always in the table output 
#' @param dropExtra Should we drop any value not specified?
#' @param defValue Default value
#' @return A named vector, the same as table but with all levels present and set to defValue if no data is present.
ltable = function(...,levels=NULL,dropExtra=TRUE,defValue=0) {
  if(is.null(levels)){
    table(...)
  }else{
    x = table(...)
    tmp = rep(defValue,length(levels))
    names(tmp)=levels
    if(dropExtra){
      o = which(names(x)%in%levels)
      tmp[names(x[o])] = x[o]
    }else{
      tmp[names(x)] = x
    }
    return(tmp)
  }
}


