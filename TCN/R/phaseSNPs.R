#' Phases SNPs using IMPUTE.
#'
#' The input is assumed to contain only SNPs heterozygous in the normal, for which we already have the core columns listed below.  IMPUTE is then run on each Region/Group combination to infer phasing. 
#' 
#' The key columns that the BAF GRanges object needs to have are:
#' \itemize{ 
#'  \item REF - The reference base.
#'  \item ALT - The alternate base.
#'  \item MUT - The number of reads covering this location in the tumour.
#'  \item MUT_ALT - The number of reads with the alternate base at this location in the tumour.
#'  \item Group - String that will be used to aggregate SNPs across different sequencing experiments (samples usually).
#'  \item Region - Will be used to aggregate SNPs at different genomic locations (typically genes).
#' }
#' And optionally:
#' \itemize{
#'  \item NORM - The number of reads covering this location in the normal.
#'  \item NORM_ALT - The number of reads with the alternate base covering this location in the normal.
#' }
#' We need to calculate:
#' \itemize{
#'  \item tau - A value giving the ratio of weights for selecting a read with the alternate and reference allele.  That is, the bias correction.
#'  \item CorrectedPhase - A series of factors (typically '0' and '1') indicating which SNPs come from the same allele.
#'  \item phaseOfMajorAllele - Used to determine which phase is the major allele (i.e., the allele with the highest CN).
#' }
#'
#' @export
#' @param bafs A GRanges object holding SNPs heterzygous in normal with columns as described above.
#' @param impute_exec The location of the IMPUTE executable.
#' @param impute_m A string giving the location of the impute m files where \%s is to be replaced by the chromosome name.
#' @param impute_h A string giving the location of the impute h files where \%s is to be replaced by the chromosome name.
#' @param impute_l A string giving the location of the impute l files where \%s is to be replaced by the chromosome name.
#' @param tgts A GRanges object with a column named regionCol which contains the GRanges spanned by the target region.  Used to provide limits to IMPUTE.  If NULL, these are set to the max/min location of SNPs in the data.
#' @param correctPhase Should we try and improve the phasing by toggling haplotype blocks that improve the variance in the estimate of both alleles CN.
#' @param tmpDir A directory we can write to where the temporary IMPUTE files will be created.
#' @param verbose Should we output diagnostic and progress messages.
#' @return bafs object with extra columns holding the phasing information
phaseSNPs = function(bafs,impute_exec,impute_m,impute_h,impute_l,tgts=NULL,correctPhase=TRUE,tmpDir='.',verbose=TRUE){
  if(length(bafs)==0){
    bafs$Phase = character()
    bafs$CorrectedPhase = character()
    bafs$Conf = numeric()
    bafs$phaseOfMajorAllele = character()
    return(bafs)
  }
  sdat = bafs$Group
  rdat = bafs$Region
  patients = unique(sdat)
  genes = unique(rdat)
  bafs$Phase = NA
  bafs$CorrectedPhase = NA
  bafs$Conf = NA
  bafs$phaseOfMajorAllele = NA
  for(patient in patients){
    for(gene in genes){
      #Get the SNPs to use
      o = which(sdat==patient & rdat==gene)
      #What chromosome is the gene on?
      chrs = unique(as.character(seqnames(bafs[o])))
      if(length(chrs)!=1){
        cat(sprintf('Region %s spans multiple chromosomes.  Leaving un-phased.\n',gene))
        bafs$Phase[o]='0'
        bafs$CorrectedPhase[o]='0'
        bafs$phaseOfMajorAllele[o] = '0'
        next
      }
      chr = chrs[1]
      #Check we have enough SNPs to bother with
      if(length(o)<2){
        bafs$Phase[o] = '0'
        bafs$CorrectedPhase[o] = '0'
        bafs$phaseOfMajorAllele[o] = '0'
        next
      }
      #Convert to the stupid format required
      gt_table = data.frame(sprintf('SNP%d',seq_along(o)),sprintf('rs%d',seq_along(o)),start(bafs[o]),bafs$REF[o],bafs$ALT[o],0,1,0)
      #Write out the genotype file in the required stupid format
      gt_file = sprintf('%s.%s.%s.gens',patient,gene,chr)
      gt_file = file.path(tmpDir,gt_file)
      write.table(gt_table,gt_file,row.names=FALSE,col.names=FALSE,quote=FALSE)
      m_file = sprintf(impute_m, chr)
      h_file = sprintf(impute_h, chr)
      l_file = sprintf(impute_l, chr)
      o_file = sprintf('%s.%s.%s.impute.phasing',patient,gene,chr)
      o_file = file.path(tmpDir,o_file)
      #Get the chromosome limits
      if(!is.null(tgts)){
        chr_max = max(end(tgts[tgts$Region == gene]))
        chr_min = min(start(tgts[tgts$Region == gene]))
      }else{
        chr_max = max(start(bafs[o]))
        chr_min = min(start(bafs[o]))
      }
      #Construct the command
      cmd = sprintf("%s -phase -m %s -h %s -l %s -g %s -int %0.0f %0.0f -o %s", impute_exec, m_file, h_file, l_file, gt_file, chr_min,chr_max,o_file)
      if(verbose){
        cat(sprintf('Phasing SNPs in %s in region %s\n',patient,gene))
      }
      system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
      #Read in the haplotypes
      haps = read.table(sprintf('%s_haps',o_file))
      haps_conf = read.table(sprintf('%s_haps_confidence',o_file))
      #Sometimes IMPUTE drops a SNP because it's low quality.  In which case we should make sure we're phasing the right ones
      oo=as.integer(gsub('SNP','',haps[,1]))
      #Work out where to put them
      #Add the phasing back into the main table
      bafs$Phase[o[oo]] = haps[,6]
      bafs$Conf[o[oo]] = haps_conf[,6]
      #Remove the files we created
      file.remove(sprintf('%s_haps',o_file))
      file.remove(sprintf('%s_haps_confidence',o_file))
      file.remove(sprintf('%s_info',o_file))
      file.remove(sprintf('%s_info_by_sample',o_file))
      file.remove(sprintf('%s_summary',o_file))
      file.remove(sprintf('%s_warnings',o_file))
      file.remove(o_file)
      file.remove(gt_file)
      #Calculate a rough estimate of the allele frequency
      est_afs = bafs$MUT_ALT[o[oo]]/bafs$MUT[o[oo]]
      #Correct it for bias if we have the bias factor
      taus = if('tau' %in% colnames(mcols(bafs))) bafs$tau[o[oo]] else rep(1,length(o))
      est_afs = est_afs/(taus*(1-est_afs)+est_afs)
      if(correctPhase){
        #Try and improve the phasing by minimising variance in each category
        phase = Rle(bafs$Phase[o[oo]])
        #Swapped phase values
        sphase = as.character(as.integer(!as.integer(as.character(runValue(phase)))))
        #Test toggling the haploblocks
        toggle = rep(FALSE,nrun(phase))
        for(j in seq_len(nrun(phase))){
          #Reset phase
          tphase = phase
          #Check if changing it reduces variance in both estimates
          runValue(tphase)[j] = sphase[j]
          toggle[j] = all(sapply(split(est_afs,as.character(tphase)),var)<sapply(split(est_afs,as.character(phase)),var))
        }
        #Actually toggle them in the copy
        tphase = phase
        runValue(tphase)[which(toggle)] = sphase[which(toggle)]
        bafs$CorrectedPhase[o[oo]] = as.character(tphase)
      }else{
        bafs$CorrectedPhase[o[oo]] = bafs$Phase[o[oo]]
      }
      #Guess which allele is the major allele
      tmp = sapply(split(est_afs,factor(bafs$CorrectedPhase[o[oo]],levels=c('0','1'))),mean)
      #If one is NA, then we have to guess based on those that aren't
      if(any(is.na(tmp))){
        if(is.na(tmp['0'])){
          maj = ifelse(tmp['1']>=0.5,'1','0')
        }else{
          maj = ifelse(tmp['0']>=0.5,'0','1')
        }
      }else{
        maj = ifelse(tmp['0']>tmp['1'],'0','1')
      }
      bafs$phaseOfMajorAllele[o[oo]] = maj
    }
  }
  return(bafs)
}
