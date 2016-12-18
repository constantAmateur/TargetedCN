#' Gets the counts of window GC content at the windows anchored on the positions given. 
#'
#' For each location referred to by the GRanges contained within reg_list, the GC content in two surrounding windows is calculated.  The first, the "forward window", start at the genomic location of interest + lstrip and is isize-lstrip-rstrip bases long.  The second, the "reverse window", ends at the genomic location of interest - lstrip and extends backwards (towards lower genomic coordinates) for isize-lstrip-rstrip bases.
#'
#' There are three different methods used to fetch this information depending on the size and length of the genomic ranges requested.  
#' When the number of genomic locations to interrogate is small and they cover only a small fraction of the genome, the sequence for each window is simply directly fetched.
#' When only a small fraction of genome is covered by a large number of independent regions, we fetch the sequence for an overlapping larger region and then extract the GC content at the locations of interest.
#' When a large fraction of the genome is required, the entire genome is fetched and GC content calculated in a sliding window.  The required windows are then extracted by sub setting this larger list.
#'
#' It is useful to specify multiple regions of interest at the same time if possible as this prevents the genome sequence having to fetched multiple times.
#'
#' After the GC has been calculated in each of the desired windows, optional aggregation is performed.  Returning either the GC count at each window, the frequency of windows with different counts for each separate GRange specified or the frequency of windows with different counts for each region in reg_list.
#'
#' @export
#' @importFrom Biostrings alphabet letterFrequencyInSlidingView
#' @importFrom GenomeInfoDb seqlevels renameSeqlevels seqlengths
#' @param reg_list A list (or GRangesList) of regions where we want to know the GC content in a window anchored at the contained locations.  For example, a list might contain a set of read locations in the first entry and a mappability mask of the genome in the second.
#' @param isize The insert size to use.
#' @param lstrip The number of bases to trim at the 5' end of the window.
#' @param rstrip The number of bases to trim at the 3' end of the window.
#' @param summarise One of 'Region', 'GRange', or 'BP'.  If 'BP', the GC count in each anchored window is returned.  If 'Region' or 'Grange', returns the frequency of window GC counts aggregated across each individual GRange entry (if 'GRange' specified) or each list entry in reg_list (if 'Region' specified).
#' @param genome The BSgenome object from which genome sequence is extracted.
#' @param min_cvg_frac In deciding how to fetch data, the number of entries in the covering GRanges objects is compared to the number in each reg_list.  If this fraction is larger than min_cvg_frac (i.e., we don't save many fetches by loading the covering regions rather then the targets directly), then the direct fetch method is used.
#' @param max_genome_frac If the requested region covers more than this fraction of the genome, then fetch the entire genome for simplicity and speed.
#' @param max_quick_regions Don't try and fetch things the quick way if we have more than this many distinct regions to fetch.
#' @return If summarise is 'Region', a list of matricies with length(reg_list) rows and columns giving the number of windows with specified GC count in that region.  Otherwise reg_list is returned with extra columns 'forw' and 'revr' added holding either the frequencies of GC windows (if summarise is 'GRange') or the individual GC counts in order (if summarise is 'BP').
get_GC_counts = function(reg_list,isize,lstrip=3,rstrip=3,summarise=c('Region','GRange','BP'),genome,min_cvg_frac=0.5,max_genome_frac=0.04,max_quick_regions=2000){
  summarise = match.arg(summarise)
  wlen = isize-rstrip-lstrip
  if(summarise=='Region'){
    #Create the output matrix
    out_for = matrix(0,nrow=length(reg_list),ncol=wlen+1)
    rownames(out_for)=names(reg_list)
    colnames(out_for)= 0:wlen
    out_rev = out_for
  }else{
    #Create place to store individual counts or summary vector if we need that
    reg_list = lapply(reg_list,function(e){
                        e$forw=vector('list',length(e))
                        e$revr=vector('list',length(e))
                        e})
  }
  #Get the chromosome lengths
  clens = seqlengths(genome)
  names(clens) = fix_chr(names(clens))
  #Get the list of chromosomes we need to look at
  chrs = unique(do.call(c,lapply(reg_list,function(e) unique(fix_chr(as.character(seqnames(e)))))))
  #Order them sensibly
  chrs = chrs[order(num_chr(chrs))]
  #The remapping of things to UCSC genome chr format 
  chrFormatMap = setNames(prefix_chr(chrs),chrs)
  #Check that we haven't asked for bad regions
  if(any(sapply(reg_list,function(e) any(start(e)<lstrip+wlen) | any(end(e) > clens[fix_chr(as.character(seqnames(e)))]-lstrip-wlen+1)))){
    stop('Location requested without any valid window.')
  }
  #Make a coverage object so we know where we need data
  cvg = reduce(GRanges(unlist(lapply(reg_list,function(e) fix_chr(as.character(seqnames(e)))),use.names=FALSE),
                       IRanges(unlist(lapply(reg_list,start),use.names=FALSE),
                               unlist(lapply(reg_list,end),use.names=FALSE))))
  #Add flank that we'll need
  start(cvg) = pmax(1,start(cvg)-lstrip-wlen+1)
  end(cvg) = pmin(clens[fix_chr(as.character(seqnames(cvg)))],end(cvg)+lstrip+wlen-1)
  #Reduce it again
  cvg = reduce(cvg)
  lcvg = length(cvg)
  #Split by (fixed up) chromosome
  cvg = split(cvg,as.character(seqnames(cvg)))
  #Check if we can do any of them super quickly
  todo = rep(TRUE,length(reg_list))
  for(i in seq_along(reg_list)){
    #What fraction of the genome is covered
    genome_frac = sum(as.numeric(width(reg_list[[i]])))/sum(as.numeric(seqlengths(genome)))
    #How long is the reduced set of regions we need to fetch compared to getting them all
    cvg_frac = lcvg/length(reg_list[[i]])
    #If we can do it quickly, do it quickly...
    #Do the quick version if we're only having to fetch 20% or fewer regions if we do it the long way and the total fraction of the genome asked for is small so we don't run into memory issues
    if(lcvg < max_quick_regions & cvg_frac > min_cvg_frac & genome_frac<max_genome_frac) {
      tmp = get_GC_counts_for_region(reg_list[[i]],isize,lstrip,rstrip,summarise=summarise!='BP',genome=genome)
      if(summarise=='Region'){
        out_for[i,] = out_for[i,] + colSums(tmp$forw)
        out_rev[i,] = out_rev[i,] + colSums(tmp$revr)
      }else if(summarise=='GRange'){
        reg_list[[i]]$forw = tmp$forw
        reg_list[[i]]$revr = tmp$revr
      }else{
        reg_list[[i]]$forw = lapply(as.list(data.frame(t(tmp$forw))),`[[`,1)
        reg_list[[i]]$revr = lapply(as.list(data.frame(t(tmp$revr))),`[[`,1)
      }
      #Now mark it as done
      todo[i]=FALSE
    }
  }
  #Process one chromosome at a time, unless we did them all quickly
  if(any(todo)){
    for(chr in chrs){
      cat(sprintf('Processing chromosome %s\n',chr))
      cat('Loading sequence.\n')
      #Work out if we should get the whole chromosome, or just the relevant subset
      chr_frac = sum(width(cvg[[chr]]))/clens[chr]
      if(chr_frac > max_genome_frac){
        genome_seq = getSeq(genome,prefix_chr(chr))
      }else{
        genome_seq = getSeq(genome,renameSeqlevels(cvg[[chr]],chrFormatMap))
      }
      #Mask anything that isn't pure ACGT in the window
      drop_bases = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='')
      drop_bases_id = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='|')
      #Get the GC content in sliding windows from the sequence
      cat('Binning GC.\n')
      if(chr_frac > max_genome_frac) {
        tmp = letterFrequencyInSlidingView(genome_seq,c('GC',drop_bases),view.width=wlen)
        cat('Re-formatting.\n')
        rm(genome_seq)
        #This chomps memory for some reason.
        gc()
        tmp = Rle(ifelse(tmp[,drop_bases_id]==0,tmp[,'G|C'],NA))
        gc()
        cat('Counting at targets.\n')
        #Get the counts at each position
        for(j in seq_along(reg_list)){
          if(!todo[j]){
            next
          }
          #Which regions do we need to get
          o = which(fix_chr(as.character(seqnames(reg_list[[j]])))==chr)
          #Don't count if this is a reverse strand position
          oof = o[as.character(strand(reg_list[[j]][o]))!='-']
          oor = o[as.character(strand(reg_list[[j]][o]))!='+']
          if(summarise=='Region'){
            #Get the windows with the right offsets
            out_for[j,] = out_for[j,] + ltable(tmp[shift(ranges(reg_list[[j]][oof]),lstrip)],levels=0:wlen)
            out_rev[j,] = out_rev[j,] + ltable(tmp[shift(ranges(reg_list[[j]][oor]),-lstrip-wlen+1)],levels=0:wlen)
          }else{
            #Get windows at each target location
            gcf = lapply(shift(ranges(reg_list[[j]][oof]),lstrip),function(e) tmp[e])
            gcr = lapply(shift(ranges(reg_list[[j]][oor]),-lstrip-wlen+1),function(e) tmp[e])
            if(summarise=='GRange'){
              #Summarise at GRange level
              reg_list[[j]]$forw[o] = lapply(gcf,ltable,levels=0:wlen)
              reg_list[[j]]$revr[o] = lapply(gcr,ltable,levels=0:wlen)
            }else{
              reg_list[[j]]$forw[o] = lapply(gcf,as.numeric)
              reg_list[[j]]$revr[o] = lapply(gcr,as.numeric)
            }
          }
        }
      }else{
        tmp = lapply(genome_seq,letterFrequencyInSlidingView,letters=c('GC',drop_bases),view.width=wlen)
        cat('Re-formatting.\n')
        gc()
        #Make it a RleList, subsetting is then easier
        tmp = RleList(lapply(tmp,function(e) Rle(ifelse(e[,drop_bases_id]==0,e[,'G|C'],NA))))
        names(tmp) = seq_along(tmp)
        gc()
        cat('Counting at targets.\n')
        for(j in seq_along(reg_list)){
          if(!todo[j]){
            next
          }
          #Which regions do we need to get
          o = which(fix_chr(as.character(seqnames(reg_list[[j]])))==chr)
          #Get the target regions.  Need to be careful about chromosome names here.
          fnomMap = setNames(fix_chr(seqlevels(reg_list[[j]])),seqlevels(reg_list[[j]]))
          ot = findOverlaps(renameSeqlevels(reg_list[[j]][o],fnomMap),cvg[[chr]])
          shits = subjectHits(ot)
          if(any(sort(queryHits(ot))!=seq(length(o)))){
            stop('This should not happen!')
          }
          #Work out the right offset for each region.
          #-start+1 gets us the window at the positions in reg_list.  Then need an extra offset to get the window we actually want.
          for_pos = shift(ranges(reg_list[[j]][o]),-start(cvg[[chr]][shits])+1+lstrip)
          rev_pos = shift(ranges(reg_list[[j]][o]),-start(cvg[[chr]][shits])+1-lstrip-wlen+1)
          #Order them sensibly to allow fast allocation
          fro = order(as.numeric(shits),start(for_pos),end(for_pos))
          rro = order(as.numeric(shits),start(rev_pos),end(rev_pos))
          #Convert them to ranges lists for easy subsetting of GC counts
          for_pos = GRanges(shits[fro],for_pos[fro])
          rev_pos = GRanges(shits[rro],rev_pos[rro])
          #Don't count if this is a reverse strand position
          oof = as.character(strand(reg_list[[j]][o[fro]]))!='-'
          oor = as.character(strand(reg_list[[j]][o[rro]]))!='+'
          if(summarise=='Region'){
            out_for[j,] = out_for[j,] + ltable(tmp[as(for_pos[oof],'RangesList')]@unlistData,levels=0:wlen)
            out_rev[j,] = out_rev[j,] + ltable(tmp[as(rev_pos[oor],'RangesList')]@unlistData,levels=0:wlen)
          }else{
            #The fast and hopefully never wrong way
            gcf = as.list(split(tmp[as(for_pos[oof],'RangesList')]@unlistData,Rle(seq_along(o[fro][oof]),width(reg_list[[j]][o[fro]][oof]))))
            gcr = as.list(split(tmp[as(rev_pos[oor],'RangesList')]@unlistData,Rle(seq_along(o[rro][oor]),width(reg_list[[j]][o[rro]][oor]))))
            if(summarise=='GRange'){
              reg_list[[j]]$forw[o[fro][oof]] = lapply(gcf,ltable,levels=0:wlen)
              reg_list[[j]]$revr[o[rro][oor]] = lapply(gcr,ltable,levels=0:wlen)
            }else{
              reg_list[[j]]$forw[o[fro][oof]] = lapply(gcf,as.numeric)
              reg_list[[j]]$revr[o[rro][oor]] = lapply(gcr,as.numeric)
            }
          }
        }
      }
    }
  }
  if(summarise=='Region'){
    return(list(forw=out_for,revr=out_rev))
  }else{
    return(reg_list)
  }
}
