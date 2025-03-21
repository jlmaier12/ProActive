% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patternMatcher.R
\name{patternMatcher}
\alias{patternMatcher}
\title{Controller function for pattern-matching}
\usage{
patternMatcher(
  pileup,
  windowSize,
  minSize,
  maxSize,
  mode,
  minContigLength,
  verbose
)
}
\arguments{
\item{pileup}{A .txt file containing mapped sequencing read coverages averaged over
100 bp windows/bins.}

\item{windowSize}{The number of basepairs to average read coverage values over.}

\item{minSize}{The minimum size (in bp) of elevation or gap patterns. Default is 10000.}

\item{maxSize}{The maximum size (in bp) of elevation or gap patterns. Default is NA
(i.e. no maximum).}

\item{mode}{Either "genome" or "metagenome".}

\item{minContigLength}{The minimum contig/chunk size (in bp) to perform pattern-matching
on. Default is 25000.}

\item{verbose}{TRUE or FALSE. Print progress messages to console. Default is TRUE.}
}
\description{
Creates the pileupSubset, representative of one contig/chunk, used as input for
each individual pattern-matching function. After the information associated with
the best match for each pattern is obtained, the pattern-match with the lowest
mean absolute difference (match-score) is used for classification.
}
\keyword{internal}
