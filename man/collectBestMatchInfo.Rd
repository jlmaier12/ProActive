% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collectBestMatchInfo.R
\name{collectBestMatchInfo}
\alias{collectBestMatchInfo}
\title{Collect information regarding the pattern-match}
\usage{
collectBestMatchInfo(pattern, pileupSubset, elevOrGap, leftRightFull)
}
\arguments{
\item{pattern}{A vector containing the values associated with the pattern-match}

\item{pileupSubset}{A subset of the pileup that pertains only to the contig/chunk
currently being assessed.}

\item{elevOrGap}{Pattern-matching on 'elevation' or 'gap' pattern.}

\item{leftRightFull}{'Left' or'Right' partial gap/elevation pattern or full
elevation/gap pattern.}
}
\description{
Make a list containing the match-score, min and max pattern-match values,
the start and stop positions of the elevated or gapped region, the elevation
ratio and the classification
}
\keyword{internal}
