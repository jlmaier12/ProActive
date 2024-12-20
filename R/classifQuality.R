#' Checks for max or minimum coverage outside match region
#'
#' @param bestMatchList
#' @param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowSize The number of base pairs to average coverage values over. For compatibility with shape-matching, window sizes can only be 200, 500, 1000 and 2000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowSize is the default.
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
classifQuality <- function (bestMatchList, pileup, windowSize, mode){
  lapply(seq_along(bestMatchList), function(i){
    bestMatchInfo <- bestMatchList[[i]]
    classification <- bestMatchInfo[[7]]
    refName <- bestMatchInfo[[8]]
    pileupSubset <- pileup[which(pileup[,1]==refName),]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    matchRegion <- c(bestMatchInfo[[4]]:bestMatchInfo[[5]])
    if(classification=="Gap") {
      covIndex <- which(pileupSubset[,2]==min(pileupSubset[,2]))[1]
      covs <- order(pileupSubset[,2], decreasing=FALSE)[1:length(matchRegion)]
      } else if (classification=="Elevation"){
        covIndex <- which(pileupSubset[,2]==max(pileupSubset[,2]))[1]
        covs <- order(pileupSubset[,2], decreasing=TRUE)[1:length(matchRegion)]
      } else {
        bestMatchList[[i]] <<- c(bestMatchInfo, NA)
        return(NULL)
      }
      if (FALSE %in% (covs %in% matchRegion) == FALSE){
        bestMatchList[[i]] <<- c(bestMatchInfo, "VeryHigh")
      }else if((covIndex %in% matchRegion) == TRUE){
        bestMatchList[[i]] <<- c(bestMatchInfo,"High")
      }else{
        bestMatchList[[i]] <<- c(bestMatchInfo, "Low")
      }
  })
    return(bestMatchList)
}
