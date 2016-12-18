#' Calculates the total negative log-likelihood from a model including both BAF and coverage data for a range of different CN states, purities and ploidies.
#'
#' For the samples specified, calculates the negative log-likelihood for the different combinations of CN values specified.
#' Essentially constructs a giant table, with each gene's CN state by allele, purity ploidy and nll for coverage and BAF models.  Will have size nrho*nploidy*(max_CN+1)*(max_CN+2)/2*nRegions.  In most cases this should not be called directly, but via fitModel.
#'
#' @export
#' @param counts The table of counts giving the number of reads with 5' fragments in regions of interest in different samples.
#' @param chi The table of bias correction factors for regions and samples of interest.
#' @param sample_correction The vector containing the sample correction factors to apply.
#' @param bafs The usual BAF GRanges object, the data from which is used to calculate likelihoods for the BAF model.
#' @param tumours The tumour samples to group together and calculate nll for.
#' @param normals The normal samples to compare against (only needed for coverage model).
#' @param disp A negative-binomial over-dispersion parameter used in the coverage model.  disp=1/psi in the usual notation.
#' @param n A range of CN states to test.  That is, nllPatient will calculate the nll for models where the CN state is set to all combinations of values given in n.
#' @param purities The range of purity values to calculate nll for.
#' @param ploidies The range of ploidy values to calculate nll for.
#' @param verbose Should we print progress and diagnostic messages?
#' @return A data.frame giving the combinations of gene,CN state, purity and ploidy and the negative log-likelihoods for each combination.
nllPatient = function(counts,chi,sample_correction,bafs,tumours,normals,disp=0,n=0:5,purities=seq(0,1,length.out=101),ploidies=seq(1,6,length.out=101),verbose=TRUE){
  samples =c(tumours,normals)
  genes = rownames(counts)
  #Exclude anything that spans chromosomes
  chrs = sapply(split((as.character(seqnames(bafs))),bafs$Region),function(e) length(unique(e)))
  genes = genes[chrs[genes]==1 | is.na(chrs[genes])]
  #Create the per gene table
  cngrid = data.frame(gene = rep(genes,each=length(n)^2),
                   nA = rep(rep(n,length(n)),length(genes)),
                   nB = rep(rep(n,each=length(n)),length(genes)),
                   nllCov = 0,
                   nllBAF = 0)
  cngrid = cngrid[cngrid$nA >= cngrid$nB,]
  #Now replicate it many times for the purity/ploidy grid
  cngrid = cngrid[rep(seq_len(nrow(cngrid)),length(purities)*length(ploidies)),]
  cngrid$purity = rep(rep(purities,length(ploidies)),each=length(genes)*length(n)*(length(n)+1)/2)
  cngrid$ploidy = rep(rep(ploidies,each=length(purities)),each=length(genes)*length(n)*(length(n)+1)/2)
  if(verbose){
    pb = txtProgressBar(style=3,min=0,max=length(purities)*length(ploidies))
    setTxtProgressBar(pb,0)
    i=0
  }
  for(purity in purities){
    for(ploidy in ploidies){
      #Get the rows we're working on
      o = which(cngrid$purity==purity & cngrid$ploidy==ploidy)
      #Now calculate the two log-likelihoods
      cngrid$nllCov[o] = nllCov(counts[genes,samples],chi[genes,samples],sample_correction[samples],tumours,normals,cngrid[o,],disp=disp)$nllCov
      cngrid$nllBAF[o] = nllBAF(bafs[bafs$Region %in% genes & bafs$Group %in% samples],cngrid[o,])$nllBAF
      if(verbose){
        i=i+1
        setTxtProgressBar(pb,i)
      }
    }
  }
  if(verbose){
    cat('\n')
  }
  return(cngrid)
}
