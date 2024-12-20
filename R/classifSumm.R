#' Summarizes best shape-matches into a table
#'
#' Summarizes the classifications made in the shape_matcher function based on which shape match achieved the lowest mean absolute difference for each pileupSubset. Outputs results in a table rather than a list.
#'
#' @param bestMatchList classifications made with shape_matcher function. classifications are stored as the first item in the bestMatchList.
#' @param pileupSubset A table containing pileupSubset names, coverages averaged over 100bp windows, and pileupSubset positions associated with mapping whole-community reads to whole-community pileupSubsets
#' @param windowSize The number of base pairs to average coverage values over. For compatibility with shape-matching, window sizes can only be 200, 500, 1000 and 2000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowSize is the default.
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
classifSumm <- function(pileup, bestMatchList, windowSize, mode){
  refName <- rep(NA, length(bestMatchList))
  confidence <- rep(NA, length(bestMatchList))
  elevRatio <- rep(NA, length(bestMatchList))
  startPos <- rep(NA, length(bestMatchList))
  endPos <- rep(NA, length(bestMatchList))
  classification <- rep(NA, length(bestMatchList))
  matchSize <- rep(NA, length(bestMatchList))
  if(length(bestMatchList)==0){
    stop("No pattern-matches detected")
  }
  lapply(seq_along(bestMatchList), function(i){
    refName[i] <<- bestMatchList[[i]][[8]]
    pileupSubset <- pileup[which(pileup[,1]==bestMatchList[[i]][[8]]),]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    confidence[i] <<- bestMatchList[[i]][[9]]
    elevRatio[i] <<- bestMatchList[[i]][[6]]
    classification[i] <<- bestMatchList[[i]][[7]]
    startPos[i] <<- pileupSubset[bestMatchList[[i]][[4]],3]
    endPos[i] <<- pileupSubset[bestMatchList[[i]][[5]],3]
    matchSize[i] <<- (length(seq(pileupSubset[bestMatchList[[i]][[4]],3], pileupSubset[bestMatchList[[i]][[5]],3], windowSize))-1) * windowSize
  })
  classifSummTable <- cbind.data.frame(refName, classification, confidence, elevRatio, startPos, endPos, matchSize)
  return(classifSummTable)
}


