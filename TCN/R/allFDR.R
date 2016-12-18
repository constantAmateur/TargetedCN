#' Performs combined inference on coverage and BAF data.
#'
#' Given the BAF based calls for difference from allelic balance, combines these with likelihood ratios for the coverage data for the same samples/genes to obtain a likelihood ratio that compares the null hypothesis of copy neutral allelic balance against the alternative of some difference.
#'
#' @param bafp Output from allBAF_FDR
#' @param patMap See help for allCoverageFDR
#' @param ... Extra arguments passed to coverageFDR
#' @return The bafp data.frame with extra columns combLR, giving the combined likelihood ratio for the patient/gene combination in a given row, and comb_pval, giving the p-value for the same.
allFDR = function(bafp,patMap,...){
  bafp$combLR = NA
  bafp$comb_pval = NA
  for(pat in unique(patMap$Patient)){
    o = which(patMap$Patient == pat)
    oo = which(bafp$Group==pat)
    if(length(oo)==0){
      next
    }
    tmp = coverageFDR(...,
                      tgt_cols = patMap$Sample[o[patMap$Type[o]=='Tumour']],
                      norm_cols = patMap$Sample[o[patMap$Type[o]!='Tumour']],
                      deviance=TRUE)
    #Now combine with the likelihood ratio from bafp
    bafp$combLR[oo] = bafp$lr[oo] + tmp[match(bafp$Region[oo],names(tmp))]
    bafp$comb_pval[oo] = pchisq(bafp$combLR[oo],df=2,lower.tail=FALSE,log.p=FALSE)
  }
  #Do FDR correction
  bafp$combFDR = p.adjust(bafp$comb_pval)
  return(bafp)
}
