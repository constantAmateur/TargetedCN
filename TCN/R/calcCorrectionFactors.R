#' Calculates a correction factor for each sample.
#'
#' Calculates a factor which can be used to normalise each sample to its ploidy, under the assumption that most of the regions specified in "null_rows" have copy number equal to the sample ploidy.
#'
#' @export
#' @param counts A matrix containing the number of fragments in each group/region, with columns giving samples and rows giving region (gene).
#' @param chi A matrix with the same shape as counts but containing bias correction factors for each combination (usually counts predicted from GC content).
#' @param null_rows The rows in counts/chi that we will assume have CN = ploidy on average.
#' @param trim The amount of data to trim when calculating the trimmed mean.  If set to NA, use median instead.
#' @return A vector named by samples giving the calculated correction factors. 
calcCorrectionFactors = function(counts,chi,null_rows,trim=NA) {
  if(missing(null_rows)){
    null_rows = seq_along(counts[,1])
  }
  cnts = counts[null_rows,]/chi[null_rows,]
  if(is.na(trim)){
    return(apply(cnts,2,median))
  }else{
    return(apply(cnts,2,mean,trim=trim))
  }
}
