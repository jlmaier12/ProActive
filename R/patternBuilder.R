#' Builds pattern-match vectors
#'
#' Builds the pattern-match (vector) associated with each contig/chunk for
#' visualization.
#'
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @keywords internal
patternBuilder <- function(pileupSubset, bestMatchList) {
  minReadCov <- bestMatchList[[2]]
  maxReadCov <- bestMatchList[[3]]
  startPos <- bestMatchList[[4]]
  endPos <- bestMatchList[[5]]
  classification <- bestMatchList[[7]]
  if (classification == "Gap") {
    if (startPos == 1) {
      pattern <- c(rep(minReadCov, endPos), rep(maxReadCov, (nrow(pileupSubset) - endPos)))
    } else if (endPos == nrow(pileupSubset)) {
      pattern <- c(rep(maxReadCov, startPos), rep(minReadCov, (nrow(pileupSubset) - startPos)))
    } else {
      matchRegion <- endPos - startPos
      pattern <- c(rep(maxReadCov, startPos), rep(minReadCov, matchRegion), rep(maxReadCov, (nrow(pileupSubset) - (matchRegion + startPos))))
    }
  } else if (classification == "NoPattern") {
    pattern <- rep(minReadCov, nrow(pileupSubset))
  } else {
    if (startPos == 1) {
      pattern <- c(rep(maxReadCov, endPos), rep(minReadCov, (nrow(pileupSubset) - endPos)))
    } else if (endPos == nrow(pileupSubset)) {
      pattern <- c(rep(minReadCov, startPos), rep(maxReadCov, (nrow(pileupSubset) - startPos)))
    } else {
      matchRegion <- endPos - startPos
      pattern <- c(rep(minReadCov, startPos), rep(maxReadCov, matchRegion), rep(minReadCov, (nrow(pileupSubset) - (matchRegion + startPos))))
    }
  }
  return(pattern)
}
