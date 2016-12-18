#' Predicts the coverage at the specified location based on GC content.
#'
#' Predicts the number of fragments overlapping each base in the desired region base on the GC content in windows anchored around each base.
#'
#' @export
#' @param tgts The GRanges object containing the regions at which we want predicted coverage.
#' @param gcBias The gcBias object for this sample.
#' @return The tgts object with a predCov column containing the predicted coverage.
predict_coverage_from_GC = function(tgts,gcBias){
  #Adjust the start/stop so that we get the whole coverage
  ltgts = tgts
  start(ltgts) = pmax(1,start(ltgts)-gcBias$isize+1)
  end(ltgts) = end(ltgts)+gcBias$isize-1
  #Reduce it, so we don't get overlapping region
  ltgts = reduce(ltgts)
  #Get the predicted read number at each location
  tmp = predict_counts_from_GC(ltgts,gcBias,'BP')
  #Construct GRanges for each "read"
  reads = GRanges(rep(seqnames(tmp),width(tmp)),
                IRanges(start=rep(start(tmp),width(tmp)) + unlist(lapply(width(tmp),seq),use.names=FALSE)-1,
                        width=gcBias$isize),
                cov = unlist(tmp$forw,use.names=FALSE)
                )
  reads = c(reads,GRanges(rep(seqnames(tmp),width(tmp)),
                IRanges(end = rep(start(tmp),width(tmp)) + unlist(lapply(width(tmp),seq),use.names=FALSE)-1,
                        width=gcBias$isize),
                cov = unlist(tmp$revr,use.names=FALSE)
                ))
 
  #Get the coverage object
  cov = coverage(reads,weight=reads$cov)
  #Finally, extract the coverage at the regions we wanted in the first instance
  tgts$predCov = as.list(cov[tgts])
  return(tgts)
}
