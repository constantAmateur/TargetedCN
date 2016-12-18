#' Creates a sunrise plot for the sample specified.
#'
#' Creates a sunrise plot comparing the total negative log-likelihood of models with different purities and ploidies.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_distiller xlab ylab labs ggtitle
#' @param purities The purities to test.
#' @param ploidies The ploidies to test.
#' @param ... Extra parameters that are passed to nllPatient which is used to calculate the nlls.
#' @return A ggplot object containing the sunrise plot.
makeSunrise = function(purities=seq(0,1,length.out=101),ploidies=seq(1,12,length.out=101),...) {
  #Get the big table of negative log-likelihoods
  if(!('n' %in% list(...))){
    n=0:(ceiling(max(ploidies)/2)+1)
    cngrid = nllPatient(purities=purities,ploidies=ploidies,n=n,...)
  }else{
    cngrid = nllPatient(purities=purities,ploidies=ploidies,...)
  }
  #For each gene/purity/ploidy combo, get the minimum nll
  df = getMinNll(cngrid)
  df = df[df$rho<1,]
  df$nll[df$nll>quantile(df$nll,0.8)]=quantile(df$nll,0.8)
  #BS to make R CMD check happy
  pl = rho = nll = NULL
  gg = ggplot(df,aes(pl,rho)) + 
        geom_raster(aes(fill=nll),interpolate=TRUE) +
        scale_fill_distiller(palette='Spectral') +
        xlab('Polidy') +
        ylab('Cell purity') +
        labs(fill='Negative log-likelihood') +
        ggtitle('Sunrise plot')
  return(gg)
}
