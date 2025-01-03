#' exampleGenomeSubsetResults
#'
#' Sample output from running `ProActive()` in genome mode
#' Report...
#'
#' @keywords internal
#' @format ## `exampleGenomeSubsetResults `
#' A list with 6 items:
#' \describe{
#'   \item{SummaryTable}{A table containing all pattern-matching classifications}
#'   \item{CleanSummaryTable}{A table containing only gap and elevation pattern-match
#' classifications (i.e. noPattern classifications removed)}
#'   \item{PatternMatches}{ A list object containing information needed to visualize the
#' pattern-matches in `plotProActiveResults()`}
#'   \item{FilteredOut}{A table containing contigs/chunks that were filtered out for being
#' too small or having too low read coverage}
#'   \item{Arguments}{A list object containing arguments used for pattern-matching (windowSize,
#' mode, chunkSize, chunkContigs)}
#'   \item{GenePredictsTable}{A list object containing arguments used for pattern-matching (windowSize,
#' mode, chunkSize, chunkContigs)}
#' }
#' @details Output from running `ProActive()` with exampleGenomePileupSubset and
#' exampleGenomegffTSVSubset in genome mode with defaults.
"exampleGenomeSubsetResults"
