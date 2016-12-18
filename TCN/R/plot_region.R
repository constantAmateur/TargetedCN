#' Plot coverage and BAF data in the target region specified.
#'
#' Plots coverage in bin_size bins, after correcting for GC, composition and position effects.
#' 
#' @export
#' @importFrom ggplot2 ggplot aes facet_grid labs geom_point geom_blank geom_smooth geom_errorbar
#' @param ltgts The set of genomic regions to plot.
#' @param sampleBAM The location of the BAM file for the Tumour sample.
#' @param bafs The BAF GRanges object.  Should only contain information for the tumour sample of interest, otherwise you'll plot the BAFs from all samples in the provided GRange.  So don't do that.
#' @param sampleGC_Bias The GCbias object for the Tumour.
#' @param sampleCorrection Correction factor for sample effects for the Tumour sample.
#' @param normalBAM The location of the BAM file for the Normal sample.
#' @param normalGC_Bias The GCbias object for the Normal.
#' @param normalCorrection Correction factor for sample effects for the Normal sample.
#' @param referenceStrandBias Factor to correct for reference strand bias.  The ratio of weights favouring the alternate strand over the reference strand.
#' @param bin_size How large should the bins it's split into be? If NA, coverage is calculated rather than counts in bins.
#' @param nbins Instead of specifying a bin size, you can give the desired numbers of bins and have the bin_size set dynamically.  If bin_size is also set to something other than NA, it specifies a floor below which the bin_size is not allowed to fall.  If set to NULL, this argument is ignored.
#' @param pos_mode What is the x-axis given to the data.  Absolute simply plots the absolute genomic locations of the data, concatenate lines all bins one after another, ordered by the target region's start positions.
#' @param coverageBiasCorrect What types of bias correction to apply to the coverage data.  A string where each letter indicates the type of correction to apply (so empty string = correct for nothing).  Options are G - GC, P - Position, S - sample.
#' @param BAF_BiasCorrect What types of bias correction to apply to the BAF data.  A string where each letter indicates the type of correction to apply.  Options are G- GC, R - Reference bias.
#' @param logCov Should the coverage data be plotted on a log2 scale?
#' @param to_plot Which data to plot and in what order.  Options are 'Tumour/Normal', 'Tumour','Normal' and, 'BAF'.
#' @param coverage_yLims The limits on the coverage y-axis.
#' @return Returns a list with the gg plot object containing the desired plot and the data.frame containing the data to be fed into it.
plot_region = function(ltgts, sampleBAM, bafs, sampleGC_Bias = NULL, sampleCorrection = 1, normalBAM = NULL, normalGC_Bias = NULL, normalCorrection = 1, referenceStrandBias = 1, bin_size = 100, nbins = 100, pos_mode = c('Concatenate', 'Absolute'), coverageBiasCorrect = 'GPS', BAF_BiasCorrect = 'GR', logCov = FALSE, to_plot = c('Tumour/Normal', 'BAF'), coverage_yLims = ifelse(c(logCov, logCov), c(-1, 1), c(0, 4))){
  pos_mode = match.arg(pos_mode)
  #Order it sensibly
  ltgts = ltgts[order(seqnames(ltgts),start(ltgts))]
  #Work out what the minimum number of possible bins is
  nbins_min = length(ltgts)
  nbins_max = sum(width(ltgts))
  #Try to get roughly the desired number of bins
  if(!is.null(nbins)){
    #Do some sanity checks
    if(nbins > nbins_max){
      warning(sprintf("Desired number of bins %d, exceeds available locations %d.  Setting bin_size to 1.",nbins,nbins_max))
      nbins = nbins_max
    }
    if(nbins < nbins_min) {
      warning(sprintf("Target region composed of %d different segments, number of requested bins, %d, is less than this.  Setting number of desired bins to %d.",nbins_min,nbins,nbins_min))
      nbins = nbins_min
    }
    #Dynamically set bin_size, but only if the dynamic size is not smaller than the user set floor
    if(is.na(bin_size) | floor(nbins_max/nbins) > bin_size){
      bin_size = floor(nbins_max/nbins)
      message(sprintf("Dynamic bin size set to %d.",bin_size))
    }
  }
  ###########
  # Coverage
  #Now get the coverage
  if(any(c('Tumour/Normal','Tumour') %in% to_plot)){
    dfTum = getNormalisedCoverage(sampleBAM,ltgts,sampleGC_Bias,sampleCorrection,bin_size=bin_size,biasCorrect=coverageBiasCorrect)
    if(logCov){
      dfTum$se = dfTum$se/dfTum$y * log2(exp(1))
      dfTum$y = log2(dfTum$y)
    }
  }
  #And that of the normal if we have one
  if(any(c('Tumour/Normal','Normal') %in% to_plot)){
    dfNorm = getNormalisedCoverage(normalBAM,ltgts,normalGC_Bias,normalCorrection,bin_size=bin_size,biasCorrect=coverageBiasCorrect)
    if(logCov){
      dfNorm$se = dfNorm$se/dfNorm$y * log2(exp(1))
      dfNorm$y = log2(dfNorm$y)
    }
  }
  #Calculate the ratio
  if('Tumour/Normal' %in% to_plot){
    #Make a ratio track
    dfRatio = dfTum
    if(logCov){
      dfRatio$y = dfTum$y - dfNorm$y
      dfRatio$se = sqrt(dfTum$se^2+dfNorm$se^2)
    }else{
      dfRatio$y = dfTum$y / dfNorm$y
      dfRatio$se = dfRatio$y * sqrt((dfTum$se/dfTum$y)^2+(dfNorm$se/dfNorm$y)^2)
    }
  }
  # END:Coverage
  ##############
  #########
  # BAFs
  if('BAF' %in% to_plot){
    #First the easy part, get the BAFs in the region, do bias correction
    data = subsetByOverlaps(bafs,ltgts)
    if(length(unique(data$Group))>1){
      warning("Multiple samples worth of BAF data provided.  You probably didn't mean to do that.")
    }
    #Map them onto the ltgts
    o = findOverlaps(data,ltgts)
    data = data[queryHits(o)]
    data$itgt = subjectHits(o)
    #Decide what bias correction to apply (if any)
    tau = rep(1,length(data))
    tau = if (grepl('G',BAF_BiasCorrect)) tau * data$GC_correction else tau
    tau = if (grepl('R',BAF_BiasCorrect)) tau * referenceStrandBias else tau
    est_afs = data$MUT_ALT/data$MUT
    est_afs = est_afs/(tau*(1-est_afs)+est_afs)
    #Make a data table
    dfBAF = data.frame(x=start(data),y=est_afs,itgt = data$itgt,altIsMajor = data$phaseOfMajorAllele==data$CorrectedPhase)
  }
  # END:BAFs
  ###########
  #Combine data frames
  df = lapply(to_plot,switch,
         Tumour=cbind(dfTum,pType='Tumour'),
         Normal=cbind(dfNorm,pType='Normal'),
         `Tumour/Normal`=cbind(dfRatio,pType='Tumour/Normal'),
         BAF=cbind(dfBAF,se=rep(0,nrow(dfBAF)),pType=rep('BAF',nrow(dfBAF)))
         )
  df = do.call(rbind,df)
  #Set the order of the plots
  df$pType = factor(df$pType,levels=to_plot)
  #Bend facet_grid into plotting things.
  if(pos_mode=='Concatenate'){
    #Convert x to position along targeted regions.  Note that some of them may over-lap, which is frustrating.  Hopefully you've picked ltgts so this is not the case.
    offsets = cumsum(c(0,width(ltgts[-length(ltgts)])))+1
    df$x = df$x - start(ltgts[df$itgt]) + offsets[df$itgt]
    xlab = 'Position in concatenated bins along ordered capture regions'
  }else if(pos_mode=='Absolute'){
    #This is fine.  Already absolute.
    xlab = 'Genomic position'
  }
  #Create a sub-data.frame which throws out the points out of the data range.
  gdf = df
  o = !(gdf$pType %in% c('Tumour/Normal','Tumour','Normal')) | (gdf$y >= coverage_yLims[1] & gdf$y <= coverage_yLims[2])
  if(nrow(gdf)-length(o)>0){
    warning(sprintf("Discarding %d bins that lie out of requested range.",nrow(gdf)-length(o)))
  }
  gdf = gdf[o,]
  #Calculate the upper/lower error-bar bounds and fix them to range limits if needed
  o = (gdf$pType %in% c('Tumour/Normal','Tumour','Normal'))
  gdf$ymin = gdf$y
  gdf$ymax = gdf$y
  gdf$ymin[o] = ifelse(gdf$y[o]-gdf$se[o] < coverage_yLims[1],
                       coverage_yLims[1],
                       gdf$y[o]-gdf$se[o])
  gdf$ymax[o] = ifelse(gdf$y[o]+gdf$se[o] > coverage_yLims[2],
                       coverage_yLims[2],
                       gdf$y[o]+gdf$se[o])
  #Filter out the NaN entries
  gdf = gdf[!is.na(gdf$y),]
  #Some BS to make R CMD check happy
  ymin = ymax = pType = x = y = altIsMajor = NULL
  #Create the canvas
  gg = ggplot() +
       facet_grid(pType ~ ., scales='free_y') +
       labs(x=xlab,y='')
  #The geom_blank statements set the y limits
  if('Tumour' %in% to_plot){
    gg = gg + 
          geom_point(data=subset(gdf,pType=='Tumour'),aes(x,y,colour=altIsMajor),size=0.6,alpha=0.6) +
          geom_errorbar(data = subset(gdf,pType='Tumour'),aes(x,y,colour=altIsMajor,ymax = ymax,ymin=ymin)) +
          geom_smooth(data = subset(gdf,pType=='Tumour'),aes(x,y),method='lm',formula = y ~ 1,se=TRUE,linetype=2,size=0.25,colour='Black') +
          geom_blank(data = data.frame(x=df$x[1],y=coverage_yLims,ptype=factor('Tumour',levels=levels(df$pType))),aes(x,y))
  }
  if('Normal' %in% to_plot){
     gg = gg + 
          geom_point(data=subset(gdf,pType=='Normal'),aes(x,y,colour=altIsMajor),size=0.6,alpha=0.6) +
          geom_errorbar(data = subset(gdf,pType='Normal'),aes(x,y,colour=altIsMajor,ymax = ymax,ymin=ymin)) +
          geom_smooth(data = subset(gdf,pType=='Normal'),aes(x,y),method='lm',formula = y ~ 1,se=TRUE,linetype=2,size=0.25,colour='Black') +
          geom_blank(data = data.frame(x=df$x[1],y=coverage_yLims,ptype=factor('Normal',levels=levels(df$pType))),aes(x,y))
  }
  if('Tumour/Normal' %in% to_plot){
    gg = gg +
      geom_point(data = subset(gdf,pType=='Tumour/Normal'),aes(x,y,colour=altIsMajor),size=.6,alpha=.6) +
      geom_errorbar(data = subset(gdf,pType='Tumour/Normal'),aes(x,y,colour=altIsMajor,ymax = ymax,ymin=ymin)) +
      geom_smooth(data = subset(gdf,pType=='Tumour/Normal'),aes(x,y),method='lm',formula = y ~ 1,se=TRUE,linetype=2,size=0.25,colour='Black') +
      geom_blank(data = data.frame(x=df$x[1],y=coverage_yLims,pType=factor('Tumour/Normal',levels=levels(df$pType))),aes(x,y))
  }
  if('BAF' %in% to_plot){
    gg = gg +
        geom_point(data = subset(gdf,pType=='BAF'),aes(x,y,colour=altIsMajor)) +
        geom_smooth(data = subset(gdf,pType=='BAF'),aes(x,y,colour=altIsMajor),method='lm',formula = y ~ 1,se=TRUE,linetype=2) +
        geom_blank(data = data.frame(x=df$x[1],y=c(0,1),pType=factor('BAF',levels=levels(df$pType))),aes(x,y))
  }
  #Return the data and the plot
  return(list(df=df,gg=gg))
}
