#' Build elevation pattern, translate it across contig and shrink the width
#'
#' Build full-size elevation pattern, run it through the elevation translator function and then shrink the pattern until it reaches the minSize
#'
#'@param windowSize
#'@param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#'@param minSize
#'@param maxSize
#'@param elevOrGap
#'@keywords internal
fullElevGap <- function(pileupSubset, windowSize, minSize, maxSize, elevOrGap){
  maxReadCov <- max(pileupSubset[,2])
  minReadCov <- min(pileupSubset[,2])
  quarterReadCov <- (abs(maxReadCov - minReadCov))/4
  maxCovSteps <- (abs(maxReadCov - (minReadCov+(quarterReadCov*3))))/2
  minCovSteps <- (abs((minReadCov+quarterReadCov) - minReadCov))/2
  maxCoverages <- seq((minReadCov+(quarterReadCov*3)), maxReadCov, maxCovSteps)
  minCoverages <- seq(minReadCov,(minReadCov+quarterReadCov), minCovSteps)
  elevLength <- ifelse((nrow(pileupSubset)-(10000/windowSize)) > (maxSize/windowSize), maxSize/windowSize, nrow(pileupSubset)-(10000/windowSize))
  nonElev <- nrow(pileupSubset)-(elevLength+(5000/windowSize))
  maxOrMin <- ifelse(elevOrGap=="Elevation", minReadCov, maxReadCov)
  maxOrMin2 <- ifelse(elevOrGap=="Elevation", maxReadCov, minReadCov)
  pattern <- c(rep(maxOrMin, 5000/windowSize), rep(maxOrMin2, elevLength), rep(maxOrMin, nonElev))
  bestMatchInfo <- collectBestMatchInfo(pattern, pileupSubset, elevOrGap, "Full")
  lapply(seq_along(maxCoverages), function(maxCov){
      bestMatchInfo <<- fullElevGapShrink(minReadCov, windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, elevOrGap)
      lapply(seq_along(minCoverages), function(minCov){
        bestMatchInfo <<- fullElevGapShrink(minCoverages[[minCov]], windowSize, maxCoverages[[maxCov]], elevLength, nonElev, bestMatchInfo, pileupSubset, minSize, elevOrGap)
        })
    })
 return(bestMatchInfo)
}
