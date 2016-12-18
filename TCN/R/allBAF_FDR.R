#' Function to calculate the FDR and p-value from the BAF data
#'
#' Given a GRanges object containing bi-allelic frequency data for heterozygous SNPs, calculates the false discovery rate and p-value for rejecting the null hypothesis that the alleles are in balance at each Region/Group (Gene/Patient) combination.
#' 
#' @param data Data must be a GRanges object with several meta-data columns.  The GRange itself gives the position of each SNP.  Two columns, Group and Region (which typically refer to a "Sample" and "Gene" respectively) give which region and group the SNP belongs to.  CorrectedPhase should give (as a character) the phase of each SNP and phaseOfMajorAllele should give which phase at that SNP belongs to the chromosome with the highest CN.  MUT_ALT should given the number of reads supporting the mutant allele for the SNP and MUT should give the total number of reads covering the SNP.  Finally, tau, should give the bias correction factor for each SNP. 
#' @param epsilon_cut The null hypothesis is that 0.5<epsilon<epsilon_cut.  Determine this by looking at the fitted values for epsilon and adjusting so we exclude the exponential at epsilon=0.5.  0.52 seems to work pretty well in most cases.
#' @param ... Passed to baf_phased_grid_ml.
#' @return A table giving the p-values, ML estimates of epsilon and FDR for each pair of "Group" and "Region" for which we have data.
allBAF_FDR = function(data,epsilon_cut=0.52,...){
  mask = as.character(interaction(data$Group,data$Region,drop=TRUE))
  fit = baf_phased_grid_ml(data,mask,left=0,right=1,verbose=TRUE,...)
  #Some of the guesses at which allele is the major allele will be wrong, fix these up
  #Correct the phase for the dodgy ones.  No need to re-do, but do need to correct epsilon to be >0.5 (you'll get the same answer if you re-run the ones that have changed, but that's just a waste of time).
  o = which(fit$par<0.5)
  data$phaseOfMajorAllele[o] = ifelse(data$phaseOfMajorAllele[o]=='1','0','1')
  fit$par[o] = 1-fit$par[o]
  #Store the results
  umask = unique(mask)
  out = data.frame(Group = gsub('^(.*?)\\.(.*?)$','\\1',umask),
                   Region = gsub('^(.*?)\\.(.*?)$','\\2',umask),
                   N = 0,
                   epsilon = NA,
                   fA = NA,
                   lr = 0,
                   pval = 1)
  #Set the number of SNPs
  cnts = table(mask)
  out$N = as.numeric(cnts[umask])
  #And the fraction that are allele A 
  tmp = sapply(split(data$phaseOfMajorAllele==data$CorrectedPhase,mask),sum)
  out$fA = as.numeric((tmp/cnts[names(tmp)])[umask])
  #Set the estimates of epsilon
  out$epsilon = as.numeric(fit$par[umask])
  #Calculate likelihood from the null model
  #To do this, trick grid_ml into just calculating likelihood
  null_eps = pmin(epsilon_cut,fit$par)
  fit_null = baf_phased_grid_ml(data,mask,left=null_eps,right=null_eps+6e-10,tol=5e-10,ngrid=4,verbose=FALSE)
  #Now calculate the likelihood ratio 
  lr = 2*(-1*fit$value+fit_null$value)
  names(lr) = names(fit$value)
  out$lr = as.numeric(lr[umask])
  out$lr[is.na(out$lr)]=0
  #And the p-value using a chi-squared distribution
  pval = pchisq(lr,df=1,lower.tail=FALSE,log.p=FALSE)
  names(pval) = names(fit$value)
  out$pval = as.numeric(pval[umask])
  out$pval[is.na(out$pval)] = 1
  #Calculate false discovery rates
  out$FDR = out$pval
  out$FDR[out$N>0] = p.adjust(out$pval[out$N>0])
  return(out)
}
