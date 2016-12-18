#' TCN: A package for estimating copy number from targeted capture data. 
#'
#' This package implements two different approaches to fitting the copy number:
#'  - A less detailed, but biologically meaningful fitting procedure, using the fitCN function, which aims to detect differences in total copy number from sample ploidy and allelic imbalance.
#'  - A full, maximum likelihood fitting procedure, where the exact, integer copy number states are inferred, along with the sample purity and ploidy of the sample.
#' Unless you have an abundance of data and a strong expectation for what the sample purity and or ploidy are in advance, it is not recommended that you attempt the full model fitting procedure.  This is because such models tend to prefer solutions with biologically unlikely copy number complexity and/or low sample purity.  The only way to prevent this is to force the fitting procedure to prefer certain solutions, thus mitigating the advantage of a completely agnostic discovery procedure, or by having enough data for the correct solution to be found.  The latter case is usually only satisfied when you have an exomes worth of data or more.
#' 
#' @docType package
#' @name TCN
#' @import edgeR
#' @import GenomicFiles
#' @import Rsamtools
#' @import GenomicRanges
#' @import BiocGenerics
#' @import methods
#' @importFrom stats dbinom median model.matrix optim p.adjust pbinom pchisq quantile rbinom setNames var
#' @importFrom utils read.table write.table setTxtProgressBar txtProgressBar
#' @importFrom IRanges IRanges subsetByOverlaps RleList
#' @importFrom S4Vectors Rle subjectHits queryHits countSubjectHits runValue runValue<- runsum nrun
NULL
