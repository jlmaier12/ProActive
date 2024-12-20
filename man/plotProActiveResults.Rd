% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotProActiveResults.R
\name{plotProActiveResults}
\alias{plotProActiveResults}
\title{Plot read coverage graphs of contigs with predicted active prophages in the whole-community read coverages}
\usage{
plotProActiveResults(
  pileup,
  ProActiveResults,
  qualityFilter = NA,
  elevFilter = NA
)
}
\arguments{
\item{pileup}{A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs}

\item{ProActiveResults}{The result from `ProActive()`}

\item{qualityFilter}{A list containing pattern information associated with all contigs containing a potential active prophage. The second, third or fourth item in the list output from ProActive()}

\item{elevFilter}{Default is NA. The user may define an elevation factor in which any plots with elevated read coverages below the elevation factor are not displayed}
}
\description{
Plot the read coverages of a contig and its associated pattern-match for each contig with a predicted active prophage.
}
\examples{
\dontrun{
plotProActiveResults(whole_commreadcovs, ProActiveResults, ProActiveResults$ExtraConfidentclassifications)
}
}
