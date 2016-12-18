#' Finds the maximum likelihood values of purity and ploidy and the CN states for each gene.
#'
#' This function attempts to determine the precise copy number state, ploidy and sample purity of each group and region (patient and gene).  This works by maximising the likelihood of a joint model combining both the coverage and BAF data.  Unless you have a large number of target regions (e.g., exome data), this estimate is likely to be biased towards solutions with unrealistically complex copy number states and low purities.
#'
#' @export
#' @param rho Starting estimate for sample purity.
#' @param pl Starting estimate for sample ploidy.
#' @param ... Parameters to be passed to nllPatient. You must not specify purities and ploidies.
#' @return A list containing the fit object and the grid of nlls for different CN states for genes at the optim values of purity and ploidy.
fitModel = function(rho=0.6,pl=2.0,...) {
  #Make a temporary function which gets the minimum negative log-likelihood for the entire model
  tmpFun = function(params,...){
    purities=params[1]
    ploidies=params[2]
    cat(sprintf('Measuring ploidy = %g, purity = %g\n',ploidies,purities))
    cngrid = nllPatient(purities=purities,ploidies=ploidies,verbose=FALSE,...)
    #Work out what the total model likelihood is
    return(getMinNll(cngrid)$nll)
  }
  #This is where all the heavy lifting happens
  fit = optim(c(rho,pl),tmpFun,method='L-BFGS-B',lower=c(0.05,1),upper=c(0.95,8),...)
  #Take the best fit values and return the table for the best purity/ploidy combo
  rho=fit$par[1]
  pl=fit$par[2]
  cngrid = nllPatient(purities=rho,ploidies=pl,verbose=FALSE,...)
  return(list(fit=fit,cngrid=cngrid))
}
