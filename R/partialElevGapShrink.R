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
#'@param pileupSubset
#'@param minSize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#'@param leftOrRight
#'@keywords internal
partialElevGapShrink <- function(minCov, windowSize, maxCov, elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, leftOrRight){
  elevRatio <- maxCov/minCov
  patternL <- c(rep(maxCov, elevLength), rep(minCov, nonElev))
  if (leftOrRight=="Left"){
      repeat {
      diff <- mean(abs(pileupSubset[,2] - patternL))
      if (diff < bestMatchInfo[[1]]){
        endPos <- (which(patternL == minCov)[1])-1
        bestMatchInfo <- list(diff, minCov, maxCov, 1, endPos, elevRatio)
      }
      if (length(which(patternL==maxCov)) <= minSize/windowSize) break
      patternL <- c(patternL[-c(1:(1000/windowSize))], rep(minCov, (1000/windowSize)))
    }
  }
  if (leftOrRight=="Right"){
    patternR <- rev(patternL)
    repeat {
      diff <- mean(abs(pileupSubset[,2] - patternR))
      if (diff < bestMatchInfo[[1]]){
        startPos <- (which(patternR == maxCov)[1])
        bestMatchInfo <- list(diff, minCov, maxCov, startPos, length(patternR), elevRatio)
      }
      if (length(which(patternR==maxCov)) < minSize/windowSize) break
      patternR <- c(rep(minCov,(1000/windowSize)),patternR[-c(((length(patternR))-((1000/windowSize)-1)):length(patternR))])
    }
  }
  #bestMatchInfo <- elevOrGapClassif(bestMatchInfo, pileupSubset, leftOrRight)
  return(bestMatchInfo)
}
