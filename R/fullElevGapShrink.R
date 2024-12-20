#' Build elevation pattern, translate it across contig and shrink the width
#'
#' Build full-size elevation pattern, run it through the elevation translator function and then shrink the pattern until it reaches the minSize
#'
#'@param minCov
#'@param windowSize
#'@param maxCov
#'@param elevLength
#'@param nonElev
#'@param bestMatchInfo
#'@param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#'@param minSize
#'@param elevOrGap
#'@keywords internal
fullElevGapShrink <- function(minCov, windowSize, maxCov, elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, elevOrGap){
  maxOrMin <- ifelse(elevOrGap=="Elevation", minCov, maxCov)
  maxOrMin2 <- ifelse(elevOrGap=="Elevation", maxCov, minCov)
  pattern <- c(rep(maxOrMin, 5000/windowSize), rep(maxOrMin2, elevLength), rep(maxOrMin, nonElev))
    repeat {
      middleRows <- which(pattern == maxOrMin2)
      if (length(middleRows) <= minSize/windowSize) break
      bestMatchInfo <- patternTranslator(pileupSubset[,2], bestMatchInfo, windowSize, pattern, elevOrGap)
      pattern <- c(pattern[-c(middleRows[2]:middleRows[(1000/windowSize)+1])], rep(maxOrMin,1000/windowSize))
    }
  return(bestMatchInfo)
}

