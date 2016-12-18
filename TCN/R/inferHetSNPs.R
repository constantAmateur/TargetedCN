#' Detect SNPs heterozygous in normal.
#'
#' Determine which SNPs from the input array are heterozygous based on the counts for reference and alternate allele in the tumour or in the normal if this information is given (i.e., if the sample has a matched normal).
#' This method works by constructing two p-values, one for getting the number of alternate reads if the SNP is homozygous for the reference and one for getting the number of reference reads if the SNP is homozygous for the alternate allele.  In an ideal world, the probability in both these cases would be 0 for any number of counts greater than 0.  However, errors cause this to be untrue.  So the p-value is calculated assuming a prescribed error rate.  Do not try and set this accurately, instead set it highly conservatively to exclude false positives.
#'
#' @export
#' @param baf The GRanges BAF object.  Should already have columns REF,ALT,MUT,MUT_ALT,NORM,NORM_ALT,Group, and Region  NORM and NORM_ALT should be NA where there is no matched normal.
#' @param seq_error_rate The assumed rate at which bases appear in error when homozygous for a particular allele.
#' @param omega The degree of reference strand bias.
#' @param fdr_cut The FDR cut-off to use to declare a SNP is not homozygous.
#' @param useNorm Should we use the counts in the normal or the mutant?  A vector of length baf or a single entry.
#' @return A logical vector of length = length(baf) indicating if a SNP is heterozygous in the normal.
inferHetSNPs = function(baf,seq_error_rate=.03,omega=1,fdr_cut=.01,useNorm=TRUE){
  if(length(useNorm)!=length(baf) & length(useNorm)!=1){
    stop('useNorm must be either length 1 or the same length as baf.')
  }
  if(length(useNorm)==1){
    useNorm = rep(useNorm,length(baf))
  }
  #Get the mutant and wt counts to use at each snp
  Ms = ifelse(useNorm,baf$NORM_ALT,baf$MUT_ALT)
  Ws = ifelse(useNorm,baf$NORM-baf$NORM_ALT,baf$MUT-baf$MUT_ALT)
  #Note we do *NOT* correct for GC bias when using this approach, as the only reads we get should be the result of sequencing errors, rather than biases in fragment selection.
  #Two null-hypotheses to test.  That the SNP is AA and RR.  This is the probability of getting M or less if AA.  The probability of getting any reads not A or R in these cases is zero, except for sequencing/mapping errors.  So we need to quantify these errors to get a p-value.
  q = pbinom(Ms,Ms+Ws,omega*(1-seq_error_rate)/((1-seq_error_rate)*omega+seq_error_rate),lower.tail=TRUE)
  #The probability of getting M or more if RR.
  p = pbinom(Ms-1,Ms+Ws,(seq_error_rate*omega)/((seq_error_rate*omega)+(1-seq_error_rate)),lower.tail=FALSE)
  #Convert to FDRs
  q = p.adjust(q)
  p = p.adjust(p)
  #Get anything that looks unlikely to be homozgyous for either allele
  return((p < fdr_cut & q < fdr_cut))
}


