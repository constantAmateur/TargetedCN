#' Calculates the minimum negative log-likelihood for each purity/ploidy value.
#'
#' Given a table with all the possible CN states, returns the values which minimise the negative log-likelihood.
#'
#' @param cngrid The data.frame returned by nllPatient.
#' @param calcWith A vector giving the different contributions to the total negative log-likelihood to sum, minimise, and plot.
#' @return A data.frame giving the minimum nll for each combination of purity and ploidy.
getMinNll = function(cngrid,calcWith=c('nllBAF','nllCov')){
  #For each gene/purity/ploidy combo, get the minimum nll
  mark = as.character(with(cngrid,interaction(gene,purity,ploidy,drop=TRUE,sep='_')))
  tmp = split(cngrid,mark)
  tmp = lapply(tmp,function(e){e[order(rowSums(e[,calcWith,drop=FALSE]),abs(e$nA+e$nB-2),e$nA-e$nB)[1],]})
  rgrid = do.call(rbind,tmp)
  #Get the total ll for each purity/ploidy combo
  mark = as.character(with(rgrid,interaction(purity,ploidy,drop=TRUE,sep='_')))
  tmp = split(rgrid,mark)
  vals = sapply(tmp,function(e){sum(rowSums(e[,calcWith,drop=FALSE]))})
  purities = as.numeric(gsub('^(.*?)_(.*)$','\\1',names(vals)))
  ploidies = as.numeric(gsub('^(.*?)_(.*)$','\\2',names(vals)))
  df = data.frame(rho=purities,pl=ploidies,nll=vals)
  return(df)
}
