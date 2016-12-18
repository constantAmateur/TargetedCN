#' Truncates or drops regions that are within buf of the start/end of chromosomes
#'
#' Removes those parts of regions that either directly refer to invalid genomic coordinates or are within buf of the start or end of a chromosome.  Mostly useful for excluding locations with bad windows for GC prediction.
#'
#' @importFrom GenomeInfoDb seqlengths
#' @param regions The GRanges object to be sanitised.
#' @param genome A BSgenome object from which sequence lengths can be extracted to test for unacceptable coordinates.
#' @param buf An integer specifying how close to the start or end of a chromosome coordinates are allowed to be.
#' @return The GRanges object given in regions with ranges trimmed or removed as needed to satisfy the requirements of the function.
sanitiseRanges = function(regions,genome,buf){
  clens = seqlengths(genome)
  names(clens) = fix_chr(names(clens))
  #Find ones that are too close to start
  o = which(start(regions) < buf)
  #Drop anything that's entirely too close
  oo = end(regions[o]) < buf
  to_drop = o[oo]
  #Edit the start of the others
  start(regions[o[!oo]]) = buf
  #Find those too close to the end
  o = which(end(regions) > clens[fix_chr(as.character(seqnames(regions)))] - buf+1)
  #Drop anything that's no good
  oo = start(regions[o]) > clens[fix_chr(as.character(seqnames(regions[o])))] - buf+1
  to_drop = c(to_drop,o[oo])
  #Fix the other
  end(regions[o[!oo]]) = clens[fix_chr(as.character(seqnames(regions[o[!oo]])))] -buf +1
  #Drop the duds and return
  if(length(to_drop)>0){
    regions[-to_drop]
  }else{
    regions
  }
}

