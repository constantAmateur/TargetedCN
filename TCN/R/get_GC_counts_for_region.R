#' Get's the GC content in windows anchored at the small number of genomic locations given
#'
#' Unlike get_GC_counts, this function is optimised to get the GC content at a small number of genomic locations.  It should probably never be directly called, instead use get_GC_counts which will call this function when appropriate.
#'
#' @importFrom Biostrings alphabet letterFrequencyInSlidingView
#' @importFrom GenomeInfoDb seqlengths
#' @param regions The GRanges object giving the positions at which GC should be calculated.
#' @param isize The insert size to use at each location (including lstrip and rstrip).
#' @param lstrip The number of bases skipped at the 5' end of the fragment.
#' @param rstrip The number of bases skipped at the 3' end of the fragment.
#' @param summarise If True, The total count of each GC value within each region is returned.  Otherwise, the raw GC values at each location are.
#' @param genome The BSgenome object from which data is extracted.
#' @return A list containing forward and reverse strand counts under as "forw" and "rev".  If summarise is False, each will contain a list of the same length as regions, named by the regions, containing the GC count based at each position in regions in turn.  If summarised, each entry is a matrix with columns giving the number of GCs and rows giving the different regions.
get_GC_counts_for_region = function(regions,isize,lstrip=3,rstrip=3,summarise=TRUE,genome) {
  wlen = isize-rstrip-lstrip
  #Build a new adjusted window object that gets the sequence we need
  windows = GRanges(prefix_chr(seqnames(regions)),IRanges(pmax(1,start(regions)-(isize-rstrip-1)),pmin(seqlengths(genome)[prefix_chr(seqnames(regions))],end(regions)+(isize-rstrip-1))))
  #Now get the sequence
  tgt_seq = getSeq(genome,windows)
  #Count GCs and bases that are bad
  drop_bases = paste0(grep('[^ACGT]',alphabet(tgt_seq),value=TRUE),collapse='')
  drop_bases_id = paste0(grep('[^ACGT]',alphabet(tgt_seq),value=TRUE),collapse='|')
  tgt_gc = lapply(tgt_seq,letterFrequencyInSlidingView,c('GC',drop_bases),view.width=wlen)
  tgt_gc = lapply(tgt_gc,function(e) ifelse(e[,drop_bases_id]==0,e[,'G|C'],NA))
  #Work out what the left and right offset from the regions we care about was (after truncating and sequence fetch)
  loff = start(regions)-start(windows)
  roff = end(windows)-end(regions)
  #Get all the forward strand GC values we can
  #Need to strip off the bits we added, but also the bits that refer to regions outside our target
  tgt_gc_for = lapply(seq_along(tgt_gc),function(i) c(tgt_gc[[i]][-(1:(loff[i]+lstrip))],rep(NA,(isize-rstrip-1)-roff[i])))
  #The +1 is because to get the last n of vector of length N is (N-n+1):N
  tgt_gc_rev = lapply(seq_along(tgt_gc),function(i) c(rep(NA,(isize-rstrip-1)-loff[i]),tgt_gc[[i]][-((length(tgt_gc[[i]])-roff[i]-lstrip+1):length(tgt_gc[[i]]))]))
  names(tgt_gc_for) = as.character(regions)
  names(tgt_gc_rev) = as.character(regions)
  if(summarise){
    #Count the members in each segment
    seg_counts_for = t(sapply(tgt_gc_for,ltable,levels=0:wlen))
    seg_counts_rev = t(sapply(tgt_gc_rev,ltable,levels=0:wlen))
    return(list(forw=seg_counts_for,revr=seg_counts_rev))
  }
  #Just return the numbers
  return(list(forw=tgt_gc_for,revr=tgt_gc_rev))
}
