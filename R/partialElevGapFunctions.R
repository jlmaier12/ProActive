#' Controller function for partial elevation/gap pattern-matching
#'
#' Builds partial elevation/gap pattern-match for patterns going off both the
#' left and right sides of the contig/chunk, shrinks the width, and collects
#' best match information
#'
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param minSize The minimum size (in bp) of elevation or gap patterns. Default is 10000.
#' @param maxSize The maximum size (in bp) of elevation or gap patterns. Default is NA
#' (i.e. no maximum).
#' @keywords internal
partialElevGap <- function(pileupSubset, windowSize, minSize, maxSize) {
  maxReadCov <- max(pileupSubset[, 2])
  minReadCov <- min(pileupSubset[, 2])
  quarterReadCov <- (abs(maxReadCov - minReadCov)) / 4
  maxCovSteps <- (abs(maxReadCov - (minReadCov + (quarterReadCov * 3)))) / 2
  minCovSteps <- (abs((minReadCov + quarterReadCov) - minReadCov)) / 2
  maxCoverages <- seq((minReadCov + (quarterReadCov * 3)), maxReadCov, maxCovSteps)
  minCoverages <- seq(minReadCov, (minReadCov + quarterReadCov), minCovSteps)
  elevLength <- ifelse((nrow(pileupSubset) - (5000 / windowSize)) > (maxSize / windowSize),
                       maxSize / windowSize, nrow(pileupSubset) - (5000 / windowSize))
  nonElev <- nrow(pileupSubset) - elevLength
  patternL <- c(rep(maxReadCov, elevLength), rep(minReadCov, nonElev))
  patternR <- rev(patternL)
  bestMatchInfoL <- collectBestMatchInfo(patternL, pileupSubset, "Elevation", "Left")
  bestMatchInfoR <- collectBestMatchInfo(patternR, pileupSubset, "Elevation", "Right")
  for (maxCov in seq_along(maxCoverages)) {
    bestMatchInfoL <- partialElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]],
                                            elevLength, nonElev, bestMatchInfoL,
                                            pileupSubset, minSize, "Left")
    bestMatchInfoR <- partialElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]],
                                            elevLength, nonElev, bestMatchInfoR,
                                            pileupSubset, minSize, "Right")
    for(minCov in seq_along(minCoverages)){
      bestMatchInfoL <- partialElevGapShrink(minCoverages[[minCov]], windowSize,
                                              maxCoverages[[maxCov]], elevLength,
                                              nonElev, bestMatchInfoL, pileupSubset,
                                              minSize, "Left")
      bestMatchInfoR <- partialElevGapShrink(minCoverages[[minCov]], windowSize,
                                              maxCoverages[[maxCov]], elevLength,
                                              nonElev, bestMatchInfoR, pileupSubset,
                                              minSize, "Right")
    }
  }
  bestMatchInfo <- elevOrGapClassif(list(bestMatchInfoL, bestMatchInfoR), pileupSubset)
  return(bestMatchInfo)
}

#' Shrink the width of partial elevation and gap patterns
#'
#'  Remove values from gapped/elevated region in the pattern-match vector
#' until it reaches the `minSize`.
#'
#' @param minCov The minimum value of the pattern-match vector.
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
#' @param leftOrRight 'Left' or 'Right' partial gap/elevation pattern.
#' @keywords internal
partialElevGapShrink <- function(minCov, windowSize, maxCov, elevLength, nonElev,
                                 bestMatchInfo, pileupSubset, minSize, leftOrRight) {
  elevRatio <- maxCov / minCov
  patternL <- c(rep(maxCov, elevLength), rep(minCov, nonElev))
  if (leftOrRight == "Left") {
    repeat {
      diff <- mean(abs(pileupSubset[, 2] - patternL))
      if (diff < bestMatchInfo[[1]]) {
        endPos <- (which(patternL == minCov)[1]) - 1
        bestMatchInfo <- list(diff, minCov, maxCov, 1, endPos, elevRatio)
      }
      if (length(which(patternL == maxCov)) <= minSize / windowSize) break
      patternL <- c(patternL[-c(1:(1000 / windowSize))], rep(minCov, (1000 / windowSize)))
    }
  }
  if (leftOrRight == "Right") {
    patternR <- rev(patternL)
    repeat {
      diff <- mean(abs(pileupSubset[, 2] - patternR))
      if (diff < bestMatchInfo[[1]]) {
        startPos <- (which(patternR == maxCov)[1])
        bestMatchInfo <- list(diff, minCov, maxCov, startPos, length(patternR), elevRatio)
      }
      if (length(which(patternR == maxCov)) < minSize / windowSize) break
      patternR <- c(rep(minCov, (1000 / windowSize)), patternR[-c(((length(patternR)) - ((1000 / windowSize) - 1)):length(patternR))])
    }
  }
  return(bestMatchInfo)
}

#' Classifies partial elevation/gap pattern-matches
#'
#' classify the contig/chunk as 'gap' if the elevated region is less than 50% of
#' the length of the contig/chunk and otherwise classify as 'elevation'.
#'
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @keywords internal
elevOrGapClassif <- function(bestMatchList, pileupSubset) {
  bestMatchInfoNew <- lapply(seq_along(bestMatchList), function(i) {
    bestMatchInfo <- bestMatchList[[i]]
    startPos <- bestMatchInfo[[4]]
    endPos <- bestMatchInfo[[5]]
    leftOrRight <- ifelse(startPos == 1, "Left", "Right")
    maxCov <- bestMatchInfo[[3]]
    minCov <- bestMatchInfo[[2]]
    if (leftOrRight == "Left") {
      if (endPos < 0.5 * (nrow(pileupSubset))) {
        bestMatchInfo[[7]] <- "Elevation"
      } else {
        bestMatchInfo[[7]] <- "Gap"
        bestMatchInfo[[4]] <- endPos
        bestMatchInfo[[5]] <- nrow(pileupSubset)
      }
    } else {
      if (startPos > 0.5 * (nrow(pileupSubset))) {
        bestMatchInfo[[7]] <- "Elevation"
      } else {
        bestMatchInfo[[7]] <- "Gap"
        bestMatchInfo[[4]] <- 1
        bestMatchInfo[[5]] <- startPos
      }
    }
    bestMatchInfo
  })
  return(bestMatchInfoNew)
}
