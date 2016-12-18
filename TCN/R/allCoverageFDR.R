#' Calculates the coverage FDR for all patients
#'
#' Given a data frame describing how samples are grouped together, calculates p-values and false discovery rates for the rejection of the null hypothesis that the copy number state in each "Region"/"Group" (Gene/Patient) pair is equal to the ploidy.
#'
#' @param patMap A data.frame with columns Patient, Sample, and Type. Type should be 'Tumour' and 'Normal'.  Patient is the group designator and Sample identifies individual samples.
#' @param ... Extra arguments passed to coverageFDR
allCoverageFDR = function(patMap,...){
  covp = list()
  for(pat in unique(patMap$Patient)){
    o = which(patMap$Patient == pat)
    tmp = topTags(coverageFDR(...,
                              tgt_cols = patMap$Sample[o[patMap$Type[o]=='Tumour']],
                              norm_cols = patMap$Sample[o[patMap$Type[o]!='Tumour']],
                              deviance=FALSE),
                  n=Inf)$table
    tmp$Group = pat
    covp[[pat]] = tmp
  }
  covp = do.call(rbind,covp)
  #Add gene column
  covp$Region = gsub('.*\\.','',rownames(covp))
  return(covp)
}


