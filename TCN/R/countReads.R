#' Count the number of fragment 5' locations falling with specified target regions.
#' 
#' This is very conservative with the reads it includes, throwing out anything that doesn't look very solid.
#' 
#' @export
#' @importFrom S4Vectors countSubjectHits
#' @param bfile The location of the BAM file to load data from.
#' @param tgts A GRanges object containing the ranges that we want to fetch counts for.
#' @return The tgts GRange with extra columns forCounts, revCounts and unknownCounts giving the number of fragments with 5' in the target region by fragment orientation.
countReads = function(bfile,tgts) {
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
    strand = pos < mpos
    GRanges(unlist(lapply(e,`[[`,'rname')),
            IRanges(ifelse(strand,pos,pos+qwidth-1),width=1),
            strand=ifelse(strand,'+','-'),
            readLen=qwidth,
            insert=unlist(lapply(e,`[[`,'isize'),use.names=FALSE)
            )
  }
  #Keep just the count in each region
  reducer = function(x,y,...) {
    #Is this the first time?
    if(class(x)!='data.frame'){
      xp = countSubjectHits(findOverlaps(x[strand(x)=='+'],tgts))
      xm = countSubjectHits(findOverlaps(x[strand(x)=='-'],tgts))
      xs = countSubjectHits(findOverlaps(x[strand(x)=='*'],tgts))
      #Make the thing to store it
      x = data.frame(row.names=as.character(tgts),forCounts = as.numeric(xp),revCounts=as.numeric(xm),unknownCounts = as.numeric(xs))
    }
    x$forCounts = x$forCounts + countSubjectHits(findOverlaps(y[strand(y)=='+'],tgts))
    x$revCounts = x$revCounts + countSubjectHits(findOverlaps(y[strand(y)=='-'],tgts))
    x$unknownCounts = x$unknownCounts + countSubjectHits(findOverlaps(y[strand(y)=='*'],tgts))
    cat(sprintf('Counted %d reads.\n',sum(x$forCounts)+sum(x$RevCounts)+sum(x$unknownCounts)))
    x
  }
  done = function(e) length(e[[1]][[1]]) == 0L
  #Actually load a subset of the reads
  cnts = reduceByYield(bf,yield,map,reducer,done)
  close(bf)
  #In case the reducer hasn't been run...
  if(class(cnts)!='data.frame'){
    cnts = reducer(cnts,GRanges())
  }
  #Merge back into
  tgts$forCounts = cnts$forCounts
  tgts$revCounts = cnts$revCounts
  tgts$unknownCounts = cnts$unknownCounts
  return(tgts)
}
