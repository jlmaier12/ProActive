#' Collect information regarding the pattern-match
#'
#' Make a list containing the match-score, min and max pattern-match values,
#' the start and stop positions of the elevated or gapped region, the elevation
#' ratio and the classification
#'
#' @param pattern A vector containing the values associated with the pattern-match
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param elevOrGap Pattern-matching on 'elevation' or 'gap' pattern.
#' @param leftRightFull 'Left' or'Right' partial gap/elevation pattern or full
#' elevation/gap pattern.
#' @keywords internal
collectBestMatchInfo <- function(pattern, pileupSubset, elevOrGap, leftRightFull) {
  diff <- mean(abs(pileupSubset[, 2] - pattern))
  elevRatio <- max(pattern) / min(pattern)
  maxOrMin <- ifelse(elevOrGap == "Elevation", min(pattern), max(pattern))
  maxOrMin2 <- ifelse(elevOrGap == "Elevation", max(pattern), min(pattern))
  endPos <- if (leftRightFull == "Left") {
    (which(pattern == maxOrMin)[1]) - 1
  } else if (leftRightFull == "Right") {
    length(pattern)
  } else {
    which(pattern == maxOrMin2)[length(which(pattern == maxOrMin2))]
  }
  startPos <- ifelse(leftRightFull == "Left", 1, (which(pattern == maxOrMin2)[1]))
  bestMatchInfo <- list(diff, min(pattern), max(pattern), startPos, endPos, elevRatio, elevOrGap)
  return(bestMatchInfo)
}
