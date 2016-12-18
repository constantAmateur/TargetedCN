#' Estimates GC bias using the method of Benjamini and Speed (NAR 2011)
#' 
#' Uses the method of Benjamini and Speed to calculate rates at which we expect reads to occur based on GC content.  Returns predicted rates at which reads to occur on the forward and reverse strand, based on their GC content in a window starting at that position of length the insert size (with lstrip and rstrip bases trimmed from the 5' and 3' end respectively).
#' 
#' As well as a BAM file to extract fragment counts from, this function also needs a GRanges object defining which regions of the genome to consider.  This should be set to those bases that have an acceptably high mappability for your data (and are in the target region if using targeted data).
#'
#' @export
#' @param bfile The path to a BAM file containing the mapped reads.
#' @param tgts A GRanges object giving the locations of the genome to consider.  
#' @param genome A BSgenome object containing the genome to consider.
#' @param lstrip The number of bases to trim from the 5' end of window in which GC is considered.
#' @param rstrip The number of bases to trim from the 3' end of window in which GC is considered.
#' @param gc_smooth Smooth the predicted rates in bins of this width.  If less than zero interpreted as fraction of insert size.
#' @param maxReads Passed to BAM file reading functions.  The maximum number of reads to load at a time to prevent memory exhaustion.
#' @return A list containing:
#' forw - A vector giving the predicted number of forward strand fragments starting at a location with the GC content in a forward window given by the vector name.
#' forw - A vector giving the predicted number of reverse strand fragments starting at a location with the GC content in a reverse window given by the vector name.
#' isize - The median insert size of the data.  The window where GC content is counted is set to this size.
#' lstrip - The value passed to this function.
#' rstrip - The value passed to this function.
#' readLen - The median read length of the data.
#' genome - The value passed to this function.
#' bfile - The value passed to this function.
#' gc_smooth - The number of GC count values over which the predicted fragment rates was smoothed.
GCnorm = function(bfile,tgts,genome,lstrip=3,rstrip=3,gc_smooth = .03,maxReads=10000000) {
  #Make a non-overlapping version of targets for 
  rtgts = reduce(tgts)
  #Filter out bad reads
  #Keep only the first mate of a read, discard anything even mildly suspicious
  flags = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
         hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
         isFirstMateRead = TRUE, isSecondMateRead = FALSE, 
         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
         isDuplicate = NA)
  bf = open(BamFile(bfile,yieldSize=maxReads))
  #Get the positions of any read in the target region
  param = ScanBamParam(flags,what=c('rname','pos','mpos','qwidth','isize'))
  yield = function(e) scanBam(e,param=param)
  map = function(e) {
      pos = unlist(lapply(e,`[[`,'pos'),use.names=FALSE)
      mpos = unlist(lapply(e,`[[`,'mpos'),use.names=FALSE)
      qwidth = unlist(lapply(e,`[[`,'qwidth'),use.names=FALSE)
      strand = pos < mpos
      tmp = GRanges(unlist(lapply(e,`[[`,'rname')),
              IRanges(ifelse(strand,pos,pos+qwidth-1),width=1),
              strand=ifelse(strand,'+','-'),
              readLen=qwidth,
              insert=unlist(lapply(e,`[[`,'isize'),use.names=FALSE)
              )
      subsetByOverlaps(tmp,rtgts)
  }
  #This is just a slightly edited REDUCEsampler
  sampleSize = maxReads
  reducer = function(x, y, ...) {
    #Haven't gotten anything useful yet...
    if(length(x)==0L & length(y)==0L){
      cat('Skipping...\n')
      return(x)
    }
    #print(x)
    yld_n = length(y)
    yld_for_n = sum(strand(y)=='+')
    yld_rev_n = sum(strand(y)=='-')
    #print(c(tot,length(x),length(y)))
    #If it's the first time, store length of x
    if(!('tot' %in% colnames(mcols(x)))) {
      x$tot = length(x)
      x$nFor = sum(strand(x)=='+')
      x$nRev = sum(strand(x)=='-')
      #If we've got options, subsample
      if(length(x) >= sampleSize){
        x = x[sample(length(x),sampleSize)]
      }
    }
    #Establish counters in y
    y$tot = numeric(length(y))+x$tot[1]
    y$nFor = numeric(length(y))+x$nFor[1]
    y$nRev = numeric(length(y))+x$nRev[1]
    #Fill up x, if it's not already filled up
    if(unique(x$tot) < sampleSize) {
      #x is established 
      #First fill up x as much as we can
      x = append(x,y)
      #If we've gone too far, trim off the fat as y
      if(length(x)>sampleSize){
        y = x[seq(sampleSize+1,length(x))]
        x = x[seq(sampleSize)]
      }else{
        #Got nothing to add then...
        y=GRanges(tot=numeric(),nFor=numeric(),nRev=numeric())
      }
    }
    #Update the counters
    x$tot = x$tot + yld_n
    x$nFor = x$nFor + yld_for_n
    x$nRev = x$nRev + yld_rev_n
    #Randomly add in some new ones from y
    keep = rbinom(1L, min(sampleSize,length(y),length(x)),min(length(x),length(y))/unique(x$tot))
    #Print
    cat(sprintf("On target reads (+,-) = %g (%g,%g)\n", unique(x$tot),unique(x$nFor),unique(x$nRev)))
    i = sample(length(x),keep)
    j = sample(length(y),keep)
    x[i] = y[j]
    x
  }
  done = function(e) length(e[[1]][[1]]) == 0L
  #Actually load a subset of the reads
  cat('Loading reads.\n')
  reads = reduceByYield(bf,yield,map,reducer,done)
  #Close the file
  close(bf)
  #If it hasn't gone through the reducer, we have all the reads
  if(length(reads)<sampleSize){
    reads = reducer(reads,GRanges())
  }
  #How many reads did we get?
  Nfor = reads$nFor[1]
  Nrev = reads$nRev[1]
  Ntot = reads$tot[1]
  #Now get a fragment size distribution
  isize = median(abs(reads$insert))
  readLen = median(reads$readLen)
  wlen = isize-lstrip-rstrip
  #Sanitise reads and regions to create a list
  reads = sanitiseRanges(reads,genome,isize-rstrip)
  rtgts = sanitiseRanges(rtgts,genome,isize-rstrip)
  reg_list = list(readSubset=reads,targetGenome=rtgts)
  #Get GC counts at locations we care about
  cnts = get_GC_counts(reg_list,isize,lstrip,rstrip,genome=genome,summarise='Region')
  #Finally calculate the GC normalisation rates
  #Smooth the result and scale appropriately
  if(gc_smooth<1){
    gc_smooth = max(1,ceiling(wlen*gc_smooth))
  }
  gc_smooth = gc_smooth + !(gc_smooth%%2)
  pad = rep(0,gc_smooth %/% 2)
  #Add a rates column to the output
  noms = c(rownames(cnts$forw),'rate')
  cnts$forw = rbind(cnts$forw,rep(0,ncol(cnts$forw)))
  rownames(cnts$forw) = noms
  #Smooth the result
  cnts$forw['rate',] = as.numeric(runsum(Rle(c(pad,cnts$forw['readSubset',],pad)),gc_smooth)) / as.numeric(runsum(Rle(c(pad,cnts$forw['targetGenome',],pad)),gc_smooth))
  #Dodgy values should be 0
  cnts$forw['rate',][!is.finite(cnts$forw['rate',])]=0
  #Scale it appropriately
  cnts$forw['rate',] = cnts$forw['rate',] * Nfor/sum(cnts$forw['rate',] * cnts$forw['targetGenome',])
  #Same for reverse
  noms = c(rownames(cnts$revr),'rate')
  cnts$revr = rbind(cnts$revr,rep(0,ncol(cnts$revr)))
  rownames(cnts$revr) = noms
  cnts$revr['rate',] = as.numeric(runsum(Rle(c(pad,cnts$revr['readSubset',],pad)),gc_smooth)) / as.numeric(runsum(Rle(c(pad,cnts$revr['targetGenome',],pad)),gc_smooth))
  cnts$revr['rate',][!is.finite(cnts$revr['rate',])]=0
  cnts$revr['rate',] = cnts$revr['rate',] * Nrev/sum(cnts$revr['rate',] * cnts$revr['targetGenome',])
  #Return everything we need to make use of this
  return(list(forw=cnts$forw,revr=cnts$revr,isize=isize,lstrip=lstrip,rstrip=rstrip,readLen=readLen,genome=genome,bfile=bfile,gc_smooth=gc_smooth))
}
