#' Predict number of fragments with 5' end in ranges from GC content
#'
#' Uses the GC bias correction given to calculate how many fragments to expect based on the GC content in windows anchored at locations in the target regions specified.
#'
#' @export
#' @param tgts A GRanges object giving the locations at which we should predict counts.
#' @param gcBias The gcBias object for this sample.
#' @param summarise Should we return the number of predicted counts per base pair ('BP'), per entry in tgts ('GRange') or in total ('Region').
#' @return If summarise=='Region' then a list giving the predicted number of forward and reverse counts.  Otherwise returns the tgts object with columns forw and revr giving the predicted counts for that particular range (if summarise is 'GRange') or at each base pair in order (if summarise is 'BP').
predict_counts_from_GC = function(tgts,gcBias,summarise=c('Region','GRange','BP')) {
  summarise = match.arg(summarise)
  tmp = get_GC_counts(list(tgts),gcBias$isize,gcBias$lstrip,gcBias$rstrip,summarise=summarise,genome=gcBias$genome)
  if(summarise=='Region'){
    return(list(forw=sum(gcBias$forw['rate',]*tmp$forw[1,]),revr=sum(gcBias$revr['rate',]*tmp$revr[1,])))
  }else if(summarise=='GRange'){
    tgts$forw = apply(tmp[[1]]$forw,1,function(e) {sum(gcBias$forw['rate',] * e)})
    tgts$revr = apply(tmp[[1]]$revr,1,function(e) {sum(gcBias$revr['rate',] * e)})
  }else{
    tgts$forw = lapply(tmp[[1]]$forw,function(e) {gcBias$forw['rate',][as.integer(e+1)]})
    tgts$revr = lapply(tmp[[1]]$revr,function(e) {gcBias$revr['rate',][as.integer(e+1)]})
  }
  return(tgts)
}
