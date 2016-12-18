#' Calculate the negative log-likelihood of different CN/purity/ploidy combinations for the coverage data. 
#'
#' This fits the glm over and over with different offsets, where those offsets are chosen to model the choice of different values of A (major copy number), B (minor copy number), purity and, ploidy.  To extract the negative log-likelihood from the model, the deviance from the model with one free parameter and with only an intercept term are compared.  What this means is that what we return is not actually the nll, but the nll minus some constant which is the nll of the null model.  For maximum likelihood estimation this should not matter, but be aware of this if you try and use it to construct a statistical test, or something else.
#' sgrid controls which CN states, purities and ploidies are tested.
#'
#' @param counts The table of counts containing the genes and samples of interest.
#' @param chi The table of bias correction factor containing the genes and samples of interest.
#' @param sample_correction The vector containing the sample correction factors to apply.
#' @param tumours The names (or column numbers in counts/chi/sample_correction) of the tumour samples to test.
#' @param normals The names (or column numbers in counts/chi/sample_correction) of the normal samples to compare against.
#' @param sgrid A data.frame containing columns gene,nA,nB,purity and, ploidy.  These specify the combinations of factors to test and calculate the negative log-likelihood for.  nA is the CN of the major allele, nB is the CN of the minor allele.
#' @param disp An over-dispersion parameter (1/psi in the usual notation) that is added to the negative binomial model.
#' @return sgrid with an extra column, nllCov, containing the negative log-likelihood estimates for each combination in sgrid.
nllCov = function(counts,chi,sample_correction,tumours,normals,sgrid,disp=0){
  samples = c(tumours,normals)
  toc = counts[,samples]
  #Calculate the table of offsets.  Remember it's a log-linear model so have to log the values
  #Only need to correct for the sample effects and let the statistical test remove any capture efficiency differences (or any other location only biases).
  too = log(chi[,samples] * outer(rep(1,nrow(toc)),sample_correction[samples]))
  #Create the common DGElist
  y = DGEList(counts=toc[,samples],group=rep(c('Tumour','Normal'),c(length(tumours),length(normals))))
  #Construct the design matrix
  groups = rep(c('Tumour','Normal'),c(length(tumours),length(normals)))
  design = model.matrix(~groups)
  design0 = design[,-2,drop=FALSE]
  #Set the over-dispersion
  y$common.dispersion = disp
  #Fit the null model
  y$offsets = too[,samples]
  fit0 = glmFit(y,design)
  fit0_null = glmFit(y,design0,prior.count=0)
  #Now adjust the offset for each entry in ngrid
  mark = as.character(with(sgrid,interaction(nA+nB,purity,ploidy,sep='_',drop=TRUE)))
  for(x in unique(mark)){
    #Work out which rows we're processing
    o = which(x==mark)
    #And the values
    CN = sgrid$nA[o[1]]+sgrid$nB[o[1]]
    rho = sgrid$purity[o[1]]
    pl = sgrid$ploidy[o[1]]
    #Calculate the bit we need to add to the offset
    Delta = log((rho*CN+(1-rho)*2)/(rho*pl+(1-rho)*2))
    if(!is.finite(Delta)){
      #A huge offset, that is equivalent to -Inf, but won't break the fit routines
      Delta=-100
    }
    #Adjust the offset
    y$offsets = too[,samples] + outer(rep(Delta,nrow(too)),samples %in% tumours)
    #This will always be effectively 0, but fit it anyway in case it is not
    fit1 = glmFit(y,design)
    #Only an intercept term, the model where the offset has to carry everything...
    fit1_null = glmFit(y,design0)
    #Now the negative log-likelihood (plus a constant) is
    sgrid$nllCov[o] = ((fit1_null$deviance - fit1$deviance)/2)[sgrid$gene[o]]
  }
  return(sgrid)
}
