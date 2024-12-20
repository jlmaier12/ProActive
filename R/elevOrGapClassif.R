#'
#'
#'
#'
#'
#'@param bestMatchInfo
#'@param pileupSubset
#'@param leftOrRight
#'@keywords internal
elevOrGapClassif <- function(bestMatchInfoList, pileupSubset){
  bestMatchInfoNew <- list()
  lapply(seq_along(bestMatchInfoList), function(i){
    bestMatchInfo <- bestMatchInfoList[[i]]
    startPos <- bestMatchInfo[[4]]
    endPos <- bestMatchInfo[[5]]
    leftOrRight <- ifelse(startPos == 1, "Left", "Right")
    maxCov <- bestMatchInfo[[3]]
    minCov <- bestMatchInfo[[2]]
    if (leftOrRight == "Left"){
      if (endPos < 0.5*(nrow(pileupSubset))){
        bestMatchInfo[[7]] <- "Elevation"
      } else {
        bestMatchInfo[[7]] <- "Gap"
        bestMatchInfo[[4]] <- endPos
        bestMatchInfo[[5]] <- nrow(pileupSubset)
      }
    } else{
      if (startPos > 0.5*(nrow(pileupSubset))){
        bestMatchInfo[[7]] <- "Elevation"
      } else {
        bestMatchInfo[[7]] <- "Gap"
        bestMatchInfo[[4]] <- 1
        bestMatchInfo[[5]] <- startPos
      }
    }
    bestMatchInfoNew[[i]] <<- bestMatchInfo
  })
  return(bestMatchInfoNew)
}
