#' Get's the coverage at the locations specified.
#'
#' Calculates the number of fragments overlapping each base specified from the paired-end data in the BAM file referenced. Should only be used for small regions, not large chunks of the genome.
#'
#' @export
#' @importFrom GenomeInfoDb seqlevels
#' @param bfile The BAM file to load data from.
#' @param tgts The GRanges at which to calculate the coverage.
#' @return The tgts object with a obsCov column containing the coverage at each location.
calculateCoverage = function(bfile,tgts) {
  #Filter out bad reads
  #Keep only the first mate of a read, discard anything even mildly suspicious
  flags = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
         hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
         isFirstMateRead = TRUE, isSecondMateRead = FALSE, 
         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
         isDuplicate = NA)
  bf = open(BamFile(bfile,yieldSize=10000000))
  #Get the positions of any read in the target region
  param = ScanBamParam(flags,what=c('rname','pos','mpos','qwidth','isize'))
  yield = function(e) scanBam(e,param=param)
  map = function(e) {
    pos = unlist(lapply(e,`[[`,'pos'),use.names=FALSE)
    mpos = unlist(lapply(e,`[[`,'mpos'),use.names=FALSE)
    qwidth = unlist(lapply(e,`[[`,'qwidth'),use.names=FALSE)
    isize = unlist(lapply(e,`[[`,'isize'),use.names=FALSE)
    strand = pos < mpos
    tmp = GRanges(unlist(lapply(e,`[[`,'rname'),use.names=FALSE),
            IRanges(ifelse(strand,pos,mpos),width=abs(isize)),
            strand = ifelse(strand,'+','-'),
            readLen = qwidth,
            insert = isize
            )
    subsetByOverlaps(tmp,tgts)
  }
  #Keep just the count in each region
  reducer = function(x,y,...) {
    c(x,y)
    ##Is this the first time?
    #if(inherits(x,'GenomicRanges')){
    #  x = coverage(x)
    #}
    ##Get the common set of sequences
    #chrs = sort(unique(c(seqlevels(y),names(x))))
    ##Add in extra chromosomes to both elements
    #seqlevels(y) = chrs
    #for(chr in chrs[!(chrs %in% names(x))]){
    #  x[[chr]] = Rle(integer())
    #}
    ##Calculate new coverage
    #y = coverage(y)
    ##Combine, ordering appropriately
    #x[chrs]+y[chrs]
  }
  done = function(e) length(e[[1]][[1]]) == 0L
  #Actually load a subset of the reads
  cat('Calculating coverage from reads.\n')
  cov = coverage(reduceByYield(bf,yield,map,reducer,done))
  close(bf)
  tgts$obsCov = as.list(cov[tgts])
  return(tgts)
}
