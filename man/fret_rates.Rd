% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fret_rates.R
\name{fret_rates}
\alias{fret_rates}
\title{Calculate lambda and fdr for each peak}
\usage{
fret_rates(peaks, segment_info, max_fdr = 0.8)
}
\arguments{
\item{peaks}{List of files with output from fret_rates_prelim}

\item{segment_info}{Maximum fdr to keep data for}

\item{max_fdr}{The permutation number to be treated as original data.
For main results, use perm=0. If you want to extract significant regions for
permuted data set perm to an integer in 1:n.perm}
}
\value{
An object that can be passed to fret_thresholds
}
\description{
Find thresholds for a range of lambda values. Calculate FDR.
}
