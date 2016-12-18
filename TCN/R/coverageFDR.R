#' Function to calculate the FDR and p-value from the coverage
#'
#' Calculates a false discovery rate for rejecting the null hypothesis that the copy number state is equal to the ploidy in each region.
#' 
#' @param counts Table of counts containing the normal and tumour samples to compare.
#' @param chi Table of bias factors containing the same rows (in the same order) as counts and containing the normal and tumour samples to compare.
#' @param sample_correction Named vector, where names are samples and values are the sample correction for composition effects.  Usually the values returned by calcCorrectionFactors.
#' @param tgt_cols The column names or numbers for the target samples.
#' @param norm_cols The column names or numbers for the normal samples.
#' @param lfc_cut The log fold-change that you consider to be equivalent to no change.  Calibrate this by looking at the logFC values coming out of the fit and excluding those likely to be normal.  To be really sure, pick out those genes with epsilon < epsilon_cut and plot their logFC (most will be CN neutral).  0.15 seems reasonable in most cases.
#' @param deviance If True, returns the deviance for the fitted model (i.e., the difference in log likelihood between the fitted and null model.
#' @return A edgeR glm fit object for each row having CN different from ploidy.
coverageFDR = function(counts,chi,sample_correction,tgt_cols,norm_cols,lfc_cut=0.15,deviance=TRUE){
  samps = c(tgt_cols,norm_cols)
  toc = counts[,samps]
  #Calculate the table of offsets.  Remember it's a log-linear model so have to log the values
  #Only need to correct for the sample effects and let the statistical test remove any capture efficiency differences (or any other location only biases).
  too = log(chi[,samps] * outer(rep(1,nrow(toc)),sample_correction[samps]))
  #Now do the edgeR test
  groups = rep(c('Tumour','Normal'),c(length(tgt_cols),length(norm_cols)))
  y = DGEList(counts=toc,group=groups)
  #Construct the design matrix
  design = model.matrix(~groups)
  #Set the over-dispersion
  y$common.dispersion = 0
  #Apply the offsets
  y$offsets = too
  #Fit the model
  fit = glmFit(y,design)
  #Can extract the deviance (likelihood ratio values) from the fit object...
  #Estimates of log(Delta) =log(mCN/mPloidy) in Tumour and Normal are in fit$coefficients
  tst = glmTreat(fit,lfc=lfc_cut)
  if(deviance){
    fit_null = glmFit(y,design[,-2,drop=FALSE])
    return(fit_null$deviance - fit$deviance)
  }
  return(tst)
}


