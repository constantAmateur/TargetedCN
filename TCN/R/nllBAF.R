#' Calculates the negative log-likelihood of the BAFs having the CN and purity values given in sgrid at the gene level.
#'
#' For each combination of gene,A,A+B,rho in sgrid, calculate the negative log-likelihood of a binomial model for the BAFs observed with these parameters.  This is nicely vectorised so, providing the memory doesn't become an issue, we can calculate all values of rho and A,A+B simultaneously.  The genes, purities and, CN states that are tested are only those that appear in sgrid.
#' We assume that all entries in bafs are from the same patient and so should be aggregated.  That is, if bafs contains SNP data from multiple samples, they will all be treated as if they came from one sample.  This is appropriate for samples which are technical or biological replicates, but obviously not for different conditions/patients.
#'
#' @param bafs The usual BAF GRanges object, the data from which is used to calculate likelihoods.
#' @param sgrid A data.frame with columns gene,nA,nB and, purity giving the combinations of gene, CN of major allele, CN of minor allele and, tumour purity to test.
#' @return sgrid with an extra column called nllBAF containing the sum of the negative log-likeihoods for the combination of parameters listed in the row summed across all SNPs in bafs for that gene.
nllBAF = function(bafs,sgrid){
  #Which genes do we need info for
  genes = unique(sgrid$gene)
  #For each SNP need to calculate one likelihood per ngrid combo (excluding ploidy)
  mark = as.character(with(sgrid,interaction(nA,nA+nB,purity,sep='_',drop=TRUE)))
  #Get the unique combinations we need to evaluate
  umark = unique(mark)
  #Work out epsilon for each
  eps = sapply(lapply(strsplit(umark,'_'),as.numeric),function(e) {(e[3]*e[1]+(1-e[3]))/(e[3]*e[2]+(1-e[3])*2)})
  #Fix up the 0,0,1 case
  eps[!is.finite(eps)]=0
  #Each SNP needs to be evaluated unique(mark) times
  nreps = length(unique(mark))
  big_altIsMajor = rep(as.character(bafs$CorrectedPhase) == as.character(bafs$phaseOfMajorAllele),nreps)
  big_A = rep(bafs$MUT_ALT,nreps)
  big_T = rep(bafs$MUT,nreps)
  big_t = rep(bafs$tau,nreps)
  big_genes = rep(bafs$Region,nreps)
  big_mark = rep(umark,each=length(bafs))
  big_eps = rep(eps,each=length(bafs))
  #Calculate the likelihoods
  big_p_maj = (big_t * big_eps) / ((big_t*big_eps) + (1-big_eps))
  big_p_min = (big_t * (1-big_eps)) / ((big_t*(1-big_eps)) + big_eps)
  #Calculate the negative log-likelihood
  nll = -1*ifelse(big_altIsMajor,
             dbinom(big_A,size=big_T,prob=big_p_maj,log=TRUE),
             dbinom(big_A,size=big_T,prob=big_p_min,log=TRUE))
  #Summarise the data by gene for each unique epsilon
  splitter = as.character(interaction(big_genes,big_mark,drop=TRUE,sep='...'))
  out = sapply(split(nll,splitter),sum)
  #Now just need to map back to sgrid
  smark = as.character(interaction(sgrid$gene,mark,sep='...',drop=TRUE))
  sgrid$nllBAF = out[smark]
  sgrid$nllBAF[is.na(sgrid$nllBAF)]=0.0
  return(sgrid)
}


