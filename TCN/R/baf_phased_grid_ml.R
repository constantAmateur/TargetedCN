#' Simultaneously calculates the maximum likelihood estimate of epsilon, the ratio of major allele CN to total CN for all groups and regions.
#'
#' This is obviously less efficient than a version based on optim on a case by case basis, but with this approach the whole problem can be vectorised, making it massively faster.  Assumes a convex shape to the likelihood function with a maximum somewhere in the range 0 < epsilon < 1.
#'
#' @param data The BAF data table.
#' @param mask A mask which tells us the rows of the data table to group together.  Usually interaction(data$Group,data$Region,drop=TRUE) or something similar.
#' @param left The starting lower-bound estimate for epsilon.  Can be the same length as data, but if of length 1 all data points get the same starting lower-bound.
#' @param right The starting upper-bound estimate for epsilon.  Can be the same length as data, but if of length 1 all data points get the same starting upper-bound.
#' @param tol Search will terminate when all values of epsilon are known to at least this accuracy.
#' @param ngrid Split the space between left and right into this many pieces each iteration.
#' @param verbose Print progress and diagnostic messages.
#' @return A list containing par, the parameter estimates and value the negative log-likelihood values of those estimates.
baf_phased_grid_ml = function(data,mask,left=0,right=1,tol=1e-6,ngrid=12,verbose=TRUE) {
    if(length(left)==1){
      left = rep(left,length(data))
    }
    if(length(right)==1){
      right = rep(right,length(data))
    }
    if(ngrid<4){
      stop('ngrid should be at least 4 to ensure convergence.')
    }
    big_altIsMajor = rep(data$CorrectedPhase==data$phaseOfMajorAllele,ngrid)
    big_A = rep(data$MUT_ALT,ngrid)
    big_T = rep(data$MUT,ngrid)
    big_t = rep(data$tau,ngrid)
    #Loop until precision is met everywhere
    i=1
    while(max(right-left)>tol){
      #Transition between left and right
      big_eps = rep(left,ngrid) +rep(seq(0,1,length.out=ngrid),each=length(data))*rep(right-left,ngrid)
      big_p_maj = (big_t * big_eps) / ((big_t*big_eps) + (1-big_eps))
      big_p_min = (big_t * (1-big_eps)) / ((big_t*(1-big_eps)) + big_eps)
      #Calculate the negative log-likelihood
      nll = -1*ifelse(big_altIsMajor,
                 dbinom(big_A,size=big_T,prob=big_p_maj,log=TRUE),
                 dbinom(big_A,size=big_T,prob=big_p_min,log=TRUE))
      #Convert to a matrix
      nll = matrix(nll,nrow=length(data),ncol=ngrid,byrow=FALSE)
      #Calculate the sum over SNPs
      nll = apply(nll,2,function(e){sapply(split(e,mask),sum,na.rm=TRUE)})
      #Work out what the maximum epsilon is for each combo
      tmp = apply(nll,1,which.min)[mask]
      #Convert to left/right values by getting matrix indicies
      left = big_eps[(pmax(1,tmp-1)-1)*length(data)+seq(length(data))]
      right = big_eps[(pmin(ngrid,tmp+1)-1)*length(data)+seq(length(data))]
      #Print some diagnostics if required
      if(verbose){
        cat(sprintf('Completed loop %d: Minimum accuracy on epsilon = %g\n',i,max(right-left)))
      }
      i=i+1
    }
    #Return the maximum likelihood parameters and values for all combos
    ml = apply(nll,1,min)[mask]
    epsilons = big_eps[(tmp-1)*length(data)+seq(length(data))]
    names(epsilons) = names(ml)
    return(list(par=epsilons,value=ml))
}
