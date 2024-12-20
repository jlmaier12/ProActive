#' Shape-matching function for active MGE in mapped read coverages
#'
#' Creates the pileupSubset, representative of one contig, that is used as input for each individual shape-building function. After the information associated with the best match for each shape is obtained, the shape with the lowest mean absolute difference (score) is chosen as the prediction for the contig being assessed.
#'
#' @param pileup A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowSize The window size used to re-average read coverage datasets
#' @param minSize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxSize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#' @param mode Either "genome" or "metagenome"
#' @param minContigLength
#'
#' @keywords internal
patternMatcher <- function (pileup, windowSize, minSize, maxSize, mode, minContigLength) {
  refNames <- unique(pileup[,1])
  bestMatchList <- list()
  filteredOutContigs <- rep(NA, length(refNames))
  reason <- rep(NA, length(refNames))
  A <- 1
  B <- 1
  C <- 1
  lapply(seq_along(refNames),function(i){
    refName <- refNames[[i]]
    if(B == floor(length(refNames)/4)){
      message("A quarter of the way done with pattern-matching")
    }
    if(B == floor(length(refNames)/2)){
      message("Half of the way done with pattern-matching")
    }
    if(B == floor((length(refNames)*3)/4)){
      message("Almost done with pattern-matching!")
    }
    B <<- B+1
    pileupSubset <- pileup[which(pileup[,1] == refName),]
     if ((nrow(pileupSubset) * 100) < minContigLength) {
     #if (pileupSubset[nrow(pileupSubset),3]< minContigLength) {
      filteredOutContigs[C] <<- refName
      reason[C] <<- "Too Short"
      C <<- C+1
      return(NULL)
    } else if (pileupSubset[(order(pileupSubset[,2], decreasing=TRUE))[50],2] <= 10) {
      filteredOutContigs[C] <<-  refName
      reason[C] <<- "Low read cov"
      C <<- C+1
      return(NULL)
    }
    pileupSubset <- changewindowSize(pileupSubset,windowSize, mode)
    noPatternBestMatch <- noPattern(pileupSubset)
    partialElevBestMatch <- partialElevGap(pileupSubset, windowSize, minSize, maxSize)
    fullElevBestMatch <- fullElevGap(pileupSubset, windowSize, minSize, maxSize, "Elevation")
    fullGapBestMatch <- fullElevGap(pileupSubset, windowSize, minSize, maxSize, "Gap")
    bestMatchSumm <- list(noPatternBestMatch,
                          partialElevBestMatch[[1]],
                          partialElevBestMatch[[2]],
                          fullElevBestMatch,
                          fullGapBestMatch)
    bestMatchScoreSumm <- c(noPatternBestMatch[[1]],
                            partialElevBestMatch[[1]][[1]],
                            partialElevBestMatch[[2]][[1]],
                            fullElevBestMatch[[1]],
                            fullGapBestMatch[[1]]) %>% as.numeric()
    bestMatch <- bestMatchSumm[[which(bestMatchScoreSumm == min(bestMatchScoreSumm))[1]]]
    bestMatchList[[A]] <<- c(bestMatch, refName)
    A <<- A+1
  })
  classifQualities <- classifQuality(bestMatchList, pileup, windowSize, mode)
  filteredOutContigsdf <- na.omit(cbind.data.frame(filteredOutContigs, reason))
  patternMatchingSumm <- list(classifQualities, filteredOutContigsdf)
  return(patternMatchingSumm)
}
