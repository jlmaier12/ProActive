#' Summarizes pattern-matching results
#'
#' Summarizes the list of pattern-matching classifications into a table.
#'
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
classifSumm <- function(pileup, bestMatchList, windowSize, mode, chunkSize) {
  if (length(bestMatchList) == 0) {
    stop("No pattern-matches detected")
  }
  refName <- vapply(seq_along(bestMatchList), function(i) {bestMatchList[[i]][[8]]}, character(1))
  elevRatio <- vapply(seq_along(bestMatchList), function(i) {bestMatchList[[i]][[6]]}, numeric(1))
  startEndPos <- lapply(seq_along(bestMatchList), function(i) {
    if(mode == "genome"){
      refName <- bestMatchList[[i]][[8]]
      chunkNumber <- as.numeric(str_extract(refName, "(?<=\\_)\\d+$")) - 1
      startPos <- (bestMatchList[[i]][[4]] * windowSize) + (chunkNumber * chunkSize)
      endPos <- (bestMatchList[[i]][[5]] * windowSize) + (chunkNumber * chunkSize)
    } else {
      startPos <- bestMatchList[[i]][[4]] * windowSize
      endPos <- bestMatchList[[i]][[5]] * windowSize
    }
    list(startPos, endPos)})
  startPos <- vapply(seq_along(startEndPos), function(i){startEndPos[[i]][[1]]}, numeric(1))
  endPos <- vapply(seq_along(startEndPos), function(i){startEndPos[[i]][[2]]}, numeric(1))
  classification <- vapply(seq_along(bestMatchList), function(i) {bestMatchList[[i]][[7]]}, character(1))
  matchSize <- vapply(seq_along(bestMatchList), function(i) {
    pileupSubset <- pileup[which(pileup[, 1] == bestMatchList[[i]][[8]]), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    (length(seq(pileupSubset[bestMatchList[[i]][[4]], 3], pileupSubset[bestMatchList[[i]][[5]], 3], windowSize)) - 1) * windowSize},
    numeric(1))
  classifSummTable <- cbind.data.frame(refName, classification, elevRatio, startPos, endPos, matchSize)
  return(classifSummTable)
}

#' Removes 'NoPattern' classifications from best match list
#'
#' Removes 'NoPattern' classifications from the list of pattern-match information
#' associated with the best pattern-matches for each contig/chunk
#'
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @keywords internal
removeNoPatterns <- function(bestMatchList) {
  newBestMatchList <- lapply(seq_along(bestMatchList), function(i) {
    bestMatchInfo <- bestMatchList[[i]]
    classification <- bestMatchInfo[[7]]
    if (classification == "NoPattern") {
      return(NULL)
    } else {
      bestMatchInfo
    }
  })
  return(return(newBestMatchList[!vapply(newBestMatchList, is.null, logical(1))]))
}
