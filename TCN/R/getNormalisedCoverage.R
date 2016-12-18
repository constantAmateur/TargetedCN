#' Calculates the normalised coverage in bins for a small set of target regions.
#'
#' Here "coverage" means "number of fragments overlapping each position" when bin_size = NA and "number of 5' fragments falling in each bin" when bin_size is a positive integer.  In both cases normalisation is done using the bias factors specified in biasCorrect.
#'
#' @param bfile The location of the BAM file to get the reads from.
#' @param ltgts A GRanges object specifying the locations to get data for.  Should be small compared to the size of the genome.
#' @param gcBias A gcBias object.
#' @param sample_correction The sample correction factor to be applied to the data.
#' @param position_correction The position corrections to apply to the data.
#' @param bin_size The size of the bins to average coverage in.
#' @param biasCorrect What type of bias correction to apply.  'G' in string = GC, 'S' in string = sample, 'P' in string = position.
#' @return A data.frame giving the normalised coverage in bins across the target region.  An altIsMajor column is set, which oscillates every time there is a change between a GRanges row.
getNormalisedCoverage = function(bfile,ltgts,gcBias,sample_correction=1,position_correction=1,bin_size=100,biasCorrect='GPS') {
  if(is.na(bin_size)){
    #We're getting true coverage, not binned counts.
    binned_reads = lapply(calculateCoverage(bfile,ltgts)$obsCov,as.numeric)
    if(grepl('G',biasCorrect)){
      binned_pred = predict_coverage_from_GC(ltgts,gcBias)$predCov
    }else{
      binned_pred = lapply(width(ltgts),function(e) Rle(values=1,lengths=e))
    }
    #Make them look similar to the binned versios
    binned_pred = lapply(binned_pred,as.numeric)
  }else{
    #Filter out bad reads
    #Keep only the first mate of a read, discard anything even mildly suspicious
    flags = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
             hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
             isFirstMateRead = TRUE, isSecondMateRead = FALSE, 
             isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
             isDuplicate = NA)
    #Get the positions of any read in the target region
    param = ScanBamParam(flags,what=c('pos','mpos'),which=ltgts)
    pos = scanBam(bfile,param=param)
    #Re-order so they're in the some order as the target regions
    o = match(as.character(ltgts),names(pos))
    pos = pos[o]
    #Add in empty lists for anything that didn't get any reads
    oo = is.na(names(pos))
    pos[oo] = vector('list',sum(oo))
    names(pos[oo]) = as.character(ltgts)[is.na(o)]
    pos[oo] = list(list(pos=integer(),mpos=integer()))
    #Get counts
    nn = sapply(pos,function(e){length(e$pos)})
    #Extract coordinates
    a = unlist(sapply(pos,'[','pos'),use.names=FALSE)
    b = unlist(sapply(pos,'[','mpos'),use.names=FALSE)
    #Infer strand
    strands = a<b
    #And 5' end of insert
    reads5p = ifelse(strands,a,a+gcBias$readLen-1)
    #Set them relative to the start of each region. +1 is for 1-indexing
    reads5p = reads5p - rep(start(ltgts),times=nn) +1
    #Re-list
    reads5p = relist(reads5p,sapply(pos,'[','pos'))
    names(reads5p) = gsub('\\.pos$','',names(reads5p))
    #Drop and bin
    binned_reads = lapply(seq_along(pos),function(i){
                     return(table(cut(reads5p[[i]],unique(c(seq(0,width(ltgts[i]),bin_size),width(ltgts[i]))),include.lowest=FALSE)))})
    names(binned_reads) = names(pos)
    #Now calculate predicted GC count in each bin
    if(grepl('G',biasCorrect)){
      gc_counts = predict_counts_from_GC(ltgts,gcBias,'BP')
      #Sum forward and reverse predictions
      tmp = relist(unlist(gc_counts$forw)+unlist(gc_counts$revr),gc_counts$forw)
      #Break into bins
      binned_pred = lapply(tmp,function(e){sapply(split(e,cut(seq_len(length(e)),unique(c(seq(0,length(e),bin_size),length(e))))),sum)})
    }else{
      #Not GC correcting, so just adjust raw counts
      binned_pred = relist(rep(1,length(unlist(binned_reads))),binned_reads)
    }
  }
  se = relist(1.96 * sqrt(unlist(binned_reads)),binned_reads)
  #Finally, calculate the normalised values
  chi = if(grepl('S',biasCorrect)) sample_correction else 1
  chi = chi* (if(grepl('P',biasCorrect)) position_correction else 1)
  ncov = relist(unlist(binned_reads)/(unlist(binned_pred)*chi),binned_reads)
  #Convert it to something plotable
  tmp = unlist(ncov)
  nbins = sapply(ncov,length)
  if(is.na(bin_size)){
    bin_mid = unlist(lapply(binned_reads,seq_along),use.names=FALSE)-1 + rep(start(ltgts)-1,times=nbins)
  }else{
    bin_left = as.numeric(gsub('^.*\\.\\((.+?),(.+?)]$','\\1',names(tmp)))
    bin_right = as.numeric(gsub('^.*\\.\\((.+?),(.+?)]$','\\2',names(tmp)))
    #Get just the middle position of the bin
    bin_mid = 0.5*(bin_left+bin_right) + rep(start(ltgts)-1,times=nbins)
  }
  #Would like to plot relative to gene model, but then we need to load the gene model...
  dfCov = data.frame(x=bin_mid,y=tmp,itgt=rep(seq_along(ncov),times=nbins))
  #Set altIsMajor to an alternating vector which switches every time we switch ltgts
  dfCov$altIsMajor = as.logical(rep(seq_along(ncov)%%2,times=nbins))
  #Calculate SE while our counts are still counts
  dfCov$se = 1.96*sqrt(unlist(binned_reads))/(unlist(binned_pred)*chi)
  return(dfCov)
}
