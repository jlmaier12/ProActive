#' Controller function for full elevation/gap pattern-matching
#'
#' Builds full elevation/gap pattern-matches, shrinks the width, and collects
#' best match information
#'
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param minSize The minimum size (in bp) of elevation or gap patterns. Default is 10000.
#' @param maxSize The maximum size (in bp) of elevation or gap patterns. Default is NA
#' (i.e. no maximum).
#' @param elevOrGap Pattern-matching on 'elevation' or 'gap' pattern.
#' @keywords internal
fullElevGap <- function(pileupSubset, windowSize, minSize, maxSize, elevOrGap) {
  maxReadCov <- max(pileupSubset[, 2])
  minReadCov <- min(pileupSubset[, 2])
  quarterReadCov <- (abs(maxReadCov - minReadCov)) / 4
  maxCovSteps <- (abs(maxReadCov - (minReadCov + (quarterReadCov * 3)))) / 2
  minCovSteps <- (abs((minReadCov + quarterReadCov) - minReadCov)) / 2
  maxCoverages <- seq((minReadCov + (quarterReadCov * 3)), maxReadCov, maxCovSteps)
  minCoverages <- seq(minReadCov, (minReadCov + quarterReadCov), minCovSteps)
  elevLength <- ifelse((nrow(pileupSubset) - (10000 / windowSize)) > (maxSize / windowSize), maxSize / windowSize, nrow(pileupSubset) - (10000 / windowSize))
  nonElev <- nrow(pileupSubset) - (elevLength + (5000 / windowSize))
  maxOrMin <- ifelse(elevOrGap == "Elevation", minReadCov, maxReadCov)
  maxOrMin2 <- ifelse(elevOrGap == "Elevation", maxReadCov, minReadCov)
  pattern <- c(rep(maxOrMin, 5000 / windowSize), rep(maxOrMin2, elevLength), rep(maxOrMin, nonElev))
  bestMatchInfo <- collectBestMatchInfo(pattern, pileupSubset, elevOrGap, "Full")
  lapply(seq_along(maxCoverages), function(maxCov) {
    bestMatchInfo <<- fullElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, elevOrGap)
    lapply(seq_along(minCoverages), function(minCov) {
      bestMatchInfo <<- fullElevGapShrink(minCoverages[[minCov]], windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, elevOrGap)
    })
  })
  return(bestMatchInfo)
}

#' Shrink the width of full elevation and gap patterns
#'
#' Remove values from gapped/elevated region in the pattern-match vector
#' until it reaches the `minSize`
#'
#' @param minCov The minimum value of the pattern-match vector.
#' @param elevLength Length of the elevated/gapped pattern-match.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param maxCov The maximum value of the pattern-match vector.
#' @param elevLength Length of the elevated/gapped pattern-match region.
#' @param nonElev Length of the non-elevated/gapped pattern-match region.
#' @param bestMatchInfo The information associated with the current best pattern-match
#' for the contig/chunk being assessed.
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param minSize The minimum size (in bp) of elevation or gap patterns. Default is 10000.
#' @param elevOrGap Pattern-matching on 'elevation' or 'gap' pattern.
#' @keywords internal
fullElevGapShrink <- function(minCov, windowSize, maxCov, elevLength, nonElev,
                              bestMatchInfo, pileupSubset, minSize, elevOrGap) {
  maxOrMin <- ifelse(elevOrGap == "Elevation", minCov, maxCov)
  maxOrMin2 <- ifelse(elevOrGap == "Elevation", maxCov, minCov)
  pattern <- c(rep(maxOrMin, 5000 / windowSize), rep(maxOrMin2, elevLength), rep(maxOrMin, nonElev))
  repeat {
    middleRows <- which(pattern == maxOrMin2)
    if (length(middleRows) <= minSize / windowSize) break
    bestMatchInfo <- patternTranslator(pileupSubset[, 2], bestMatchInfo, windowSize, pattern, elevOrGap)
    pattern <- c(pattern[-c(middleRows[2]:middleRows[(1000 / windowSize) + 1])], rep(maxOrMin, 1000 / windowSize))
  }
  return(bestMatchInfo)
}

#' Full elevation/gap pattern translator
#'
#' Translates full elevation/gap patterns across contigs/chunks 1000bp at a time.
#' Translation stops when the elevation pattern is 5000bp from the end of the contig/chunk.
#'
#' @param contigCov The read coverages that pertain to the pileupSubset
#' @param bestMatchInfo The information associated with the current best pattern-match
#' for the contig/chunk being assessed.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param pattern A vector containing the values associated with the pattern-match
#' @param elevOrGap Pattern-matching on 'elevation' or 'gap' pattern.
#' @keywords internal
patternTranslator <- function(contigCov, bestMatchInfo, windowSize, pattern, elevOrGap) {
  minPattern <- min(pattern)
  maxPattern <- max(pattern)
  minOrMax <- ifelse(elevOrGap == "Elevation", min(pattern), max(pattern))
  maxOrMin2 <- ifelse(elevOrGap == "Elevation", max(pattern), min(pattern))
  repeat {
    pattern <- c(rep(minOrMax, (1000 / windowSize)), pattern[-c((length(pattern) - ((1000 / windowSize) - 1)):length(pattern))])
    if (elevOrGap == "Elevation") {
      if (pattern[length(pattern) - (5000 / windowSize)] > minPattern) break
    } else {
      if (pattern[length(pattern) - (5000 / windowSize)] < maxPattern) break
    }
    diff <- mean(abs(contigCov - pattern))
    if (diff < bestMatchInfo[[1]]) {
      elevRatio <- max(pattern) / min(pattern)
      startPos <- which(pattern == maxOrMin2)[1]
      endPos <- which(pattern == maxOrMin2)[length(which(pattern == maxOrMin2))]
      bestMatchInfo <- list(diff, minPattern, maxPattern, startPos, endPos, elevRatio, elevOrGap)
    }
  }
  return(bestMatchInfo)
}

