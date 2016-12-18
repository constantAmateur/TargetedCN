#' Calculates GC bias correction factors for each SNP
#'
#' Uses GC bias curve to calculate the correction that should be applied in order to correct for the preference for the alternate/reference allele when the GC content changes between them.  That is, if the reference and alternate allele of a SNP have different GC content (e.g., G -> T), calculates a per-SNP correction factor to account for the change in predicted coverage for the alternate allele owing to the different GC content of the two alleles.
#'
#' @export
#' @param bafs GRanges object containing data on the SNPs for the sample of interest.
#' @param gcBias The output of GCnorm for the sample of interest.
#' @return The bafs object with an extra GC_correction column.
calculateBAF_GC_Correction = function(bafs,gcBias) {
  if(length(bafs)==0){
    bafs$GC_correction = numeric()
    return(bafs)
  }
  bafs$GC_correction = 1
  isize = gcBias$isize
  read_len = gcBias$readLen
  lstrip = gcBias$lstrip
  rstrip = gcBias$rstrip
  #Only need to correct those with a GC difference
  gc_flip = xor((bafs$REF=='G' | bafs$REF=='C'),(bafs$ALT=='G' | bafs$ALT=='C'))
  o = which(gc_flip)
  #To get coverage at the SNP, we need a read to fall either no more than read_len-1 bases before if a forward strand read or no more than read_len-1 after if it's a reverse read.  If the 5' end of the read is insert size away in the right direction we can also get coverage, but this won't change when the SNP's GC status changes so we can ignore this
  regions = GRanges((seqnames(bafs)),IRanges(start(bafs)-read_len+1,width=2*read_len-1))[o]
  #Get the GC counts at each location
  counts = get_GC_counts(list(regions),gcBias$isize,gcBias$lstrip,gcBias$rstrip,summarise='BP',genome=gcBias$genome)[[1]]
  #Get the predicted counts in the forward direction
  #For reference
  cf = sapply(counts$forw,function(e){sum(gcBias$forw['rate',][e+1])})
  #And alternate
  mcf = sapply(seq_along(o),function(e){
               sum(gcBias$forw['rate',][counts$forw[[e]]+1+ifelse(bafs[o[e]]$ALT %in% c('G','C'),1,-1)])
                   })
  #Repeat for reverse
  cr = sapply(counts$revr,function(e){sum(gcBias$revr['rate',][e+1])})
  mcr = sapply(seq_along(o),function(e){
               sum(gcBias$revr['rate',][counts$revr[[e]]+1+ifelse(bafs[o[e]]$ALT %in% c('G','C'),1,-1)])
                   })
  #counts = predict_counts_from_GC(regions,gcBias,'GRange')
  #cf = counts$forw+counts$revr
  ##Predicted coverage for reference and alt on forward strand
  #cf = sapply(counts$counts_for,function(e){sum(gcBias$forw['rate',][as.character(e)])})
  #mcf = sapply(seq_along(o),function(e) {sum(gcBias$forw['rate',][as.character(counts$counts_for[[e]]+ifelse(bafs[o[e]]$ALT %in% c('G','C'),1,-1))])})
  ##Same for reverse
  #cr = sapply(counts$counts_rev,function(e){sum(gcBias$revr['rate',][as.character(e)])})
  #mcr = sapply(seq_along(o),function(e) {sum(gcBias$revr['rate',][as.character(counts$counts_rev[[e]]+ifelse(bafs[o[e]]$ALT %in% c('G','C'),1,-1))])})
  #Correction factor is predicted alternate reference count/ predicted normal count
  correction = (mcr+mcf)/(cr+cf)
  bafs$GC_correction[o] = correction
  return(bafs)
}
