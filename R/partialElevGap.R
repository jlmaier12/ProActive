
#' Prophage off-left pattern-matcher for whole-community read coverage
#'
#' Build and translate and prophage pattern going off the left side of the graph across a contig. Stop translating when the pattern left on the contig is 10,000bp. Translate the pattern 2000bp at a time.
#' Elevation pattern off left or right side
#'
#' @param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage datasets
#' @param minSize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxSize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#' @param elevOrGap
#' @param leftRight
#' @keywords internal
partialElevGap <- function (pileupSubset, windowSize, minSize, maxSize) {
  maxReadCov <- max(pileupSubset[,2])
  minReadCov <- min(pileupSubset[,2])
  quarterReadCov <- (abs(maxReadCov - minReadCov))/4
  maxCovSteps <- (abs(maxReadCov - (minReadCov+(quarterReadCov*3))))/2
  minCovSteps <- (abs((minReadCov+quarterReadCov) - minReadCov))/2
  maxCoverages <- seq((minReadCov+(quarterReadCov*3)), maxReadCov, maxCovSteps)
  minCoverages <- seq(minReadCov,(minReadCov+quarterReadCov), minCovSteps)
  elevLength <- ifelse ((nrow(pileupSubset)-(5000/windowSize))>(maxSize/windowSize),maxSize/windowSize,nrow(pileupSubset)-(5000/windowSize))
  nonElev <- nrow(pileupSubset)-elevLength
  patternL <- c(rep(maxReadCov, elevLength), rep(minReadCov, nonElev))
  patternR <- rev(patternL)
  bestMatchInfoL <- collectBestMatchInfo(patternL, pileupSubset, "Elevation", "Left")
  bestMatchInfoR <- collectBestMatchInfo(patternR, pileupSubset, "Elevation", "Right")
  lapply(seq_along(maxCoverages), function(maxCov){
    bestMatchInfoL <<- partialElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfoL, pileupSubset, minSize, "Left")
    bestMatchInfoR <<- partialElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfoR, pileupSubset, minSize, "Right")
    lapply(seq_along(minCoverages), function(minCov){
      bestMatchInfoL <<- partialElevGapShrink(minCoverages[[minCov]], windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfoL, pileupSubset, minSize, "Left")
      bestMatchInfoR <<- partialElevGapShrink(minCoverages[[minCov]], windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfoR, pileupSubset, minSize, "Right")
      })
  })
  bestMatchInfo <- elevOrGapClassif(list(bestMatchInfoL, bestMatchInfoR), pileupSubset)
  return(bestMatchInfo)
}


