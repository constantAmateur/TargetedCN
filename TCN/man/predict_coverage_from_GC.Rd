% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_coverage_from_GC.R
\name{predict_coverage_from_GC}
\alias{predict_coverage_from_GC}
\title{Predicts the coverage at the specified location based on GC content.}
\usage{
predict_coverage_from_GC(tgts, gcBias)
}
\arguments{
\item{tgts}{The GRanges object containing the regions at which we want predicted coverage.}

\item{gcBias}{The gcBias object for this sample.}
}
\value{
The tgts object with a predCov column containing the predicted coverage.
}
\description{
Predicts the number of fragments overlapping each base in the desired region base on the GC content in windows anchored around each base.
}

