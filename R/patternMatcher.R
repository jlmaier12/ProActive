#' Controller function for pattern-matching
#'
#' Creates the pileupSubset, representative of one contig/chunk, used as input for
#' each individual pattern-matching function. After the information associated with
#' the best match for each pattern is obtained, the pattern-match with the lowest
#' mean absolute difference (match-score) is used for classification.
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param windowSize The number of basepairs to average read coverage values over.
#' @param minSize The minimum size (in bp) of elevation or gap patterns. Default is 10000.
#' @param maxSize The maximum size (in bp) of elevation or gap patterns. Default is NA
#' (i.e. no maximum).
#' @param mode Either "genome" or "metagenome".
#' @param minContigLength The minimum contig/chunk size (in bp) to perform pattern-matching
#' on. Default is 25000.
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @importFrom stats na.omit
#' @keywords internal
patternMatcher <- function(pileup, windowSize, minSize, maxSize, mode, minContigLength, verbose) {
  refNames <- unique(pileup[, 1])
  bestMatchList <- vector(mode='list', length=length(refNames))
  filteredOutContigs <- rep(NA, length(refNames))
  reason <- rep(NA, length(refNames))
  A <- 1
  B <- 1
  C <- 1
  for (i in seq_along(refNames)) {
    refName <- refNames[[i]]
    if(verbose){
    if (B == floor(length(refNames) / 4)) {
      message("A quarter of the way done with pattern-matching")
    }
    if (B == floor(length(refNames) / 2)) {
      message("Half of the way done with pattern-matching")
    }
    if (B == floor((length(refNames) * 3) / 4)) {
      message("Almost done with pattern-matching!")
    }
    B <- B + 1
    }
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    if ((nrow(pileupSubset) * 100) < minContigLength) {
      filteredOutContigs[C] <- refName
      reason[C] <- "Too Short"
      C <- C + 1
      next
    } else if (pileupSubset[(order(pileupSubset[, 2], decreasing = TRUE))[50], 2] <= 10) {
      filteredOutContigs[C] <- refName
      reason[C] <- "Low read cov"
      C <- C + 1
      next
    }
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    noPatternBestMatch <- noPattern(pileupSubset)
    partialElevBestMatch <- partialElevGap(pileupSubset, windowSize, minSize, maxSize)
    fullElevBestMatch <- fullElevGap(pileupSubset, windowSize, minSize, maxSize, "Elevation")
    fullGapBestMatch <- fullElevGap(pileupSubset, windowSize, minSize, maxSize, "Gap")
    bestMatchSumm <- list(
      noPatternBestMatch,
      partialElevBestMatch[[1]],
      partialElevBestMatch[[2]],
      fullElevBestMatch,
      fullGapBestMatch
    )
    bestMatchScoreSumm <- c(
      noPatternBestMatch[[1]],
      partialElevBestMatch[[1]][[1]],
      partialElevBestMatch[[2]][[1]],
      fullElevBestMatch[[1]],
      fullGapBestMatch[[1]]
    )
    bestMatch <- bestMatchSumm[[which(bestMatchScoreSumm == min(bestMatchScoreSumm))[1]]]
    bestMatchList[[A]] <- c(bestMatch, refName)
    A <- A + 1
  }
  bestMatchList <- (bestMatchList[!vapply(bestMatchList, is.null, logical(1))])
  filteredOutContigsdf <- na.omit(cbind.data.frame(filteredOutContigs, reason))
  patternMatchingSumm <- list(bestMatchList, filteredOutContigsdf)
  return(patternMatchingSumm)
}
