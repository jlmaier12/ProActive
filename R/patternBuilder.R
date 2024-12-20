#' pattern-builder for patterns detected in whole-community read coverages
#'
#' Builds the pattern (vector) associated with the 'best pattern-match' information associated with each contig so it can be plotted.
#'
#' @param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param classifList A list containing pattern information associated with all contigs classified by ProActive
#' @param i The index position associated with the current contig's best pattern-match information
#'
#' @keywords internal
patternBuilder <- function(pileupSubset, classifList){
  minReadCov <- classifList[[2]]
  maxReadCov <- classifList[[3]]
  startPos <- classifList[[4]]
  endPos <- classifList[[5]]
  classification <- classifList[[7]]
  if(classification == "Gap") {
    if (startPos==1) {
      pattern <- c(rep(minReadCov,endPos), rep(maxReadCov, (nrow(pileupSubset)-endPos)))
    } else if (endPos == nrow(pileupSubset)){
      pattern <- c(rep(maxReadCov, startPos), rep(minReadCov, (nrow(pileupSubset)-startPos)))
    } else{
      matchRegion <- endPos-startPos
      pattern <- c(rep(maxReadCov, startPos), rep(minReadCov, matchRegion), rep(maxReadCov, (nrow(pileupSubset)-(matchRegion+startPos))))
    }
  } else if (classification == "NoPattern"){
      pattern <- rep(minReadCov, nrow(pileupSubset))
  } else {
  if (startPos==1) {
    pattern <- c(rep(maxReadCov,endPos), rep(minReadCov, (nrow(pileupSubset)-endPos)))
  } else if (endPos == nrow(pileupSubset)){
    pattern <- c(rep(minReadCov, startPos), rep(maxReadCov, (nrow(pileupSubset)-startPos)))
  } else{
    matchRegion <- endPos-startPos
    pattern <- c(rep(minReadCov, startPos), rep(maxReadCov, matchRegion), rep(minReadCov, (nrow(pileupSubset)-(matchRegion+startPos))))
  }
  }
  return(pattern)
}
