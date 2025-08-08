#' Change the pileup window size
#'
#' Re-averages windows of pileup files with 100bp windows to reduce pileup size.
#'
#' @param pileupSubset A subset of the pileup that pertains only to the contig/chunk
#'  currently being assessed.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param chunkContigs TRUE or FALSE, If TRUE and `mode`="metagenome", contigs longer
#' than the `chunkSize` will be 'chunked' into smaller subsets and pattern-matching
#' will be performed on each subset. Default is FALSE.
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
changewindowSize <- function(pileupSubset, windowSize, chunkContigs, mode) {
  coverage <- vector()
  X <- 0
  Y <- windowSize / 100
  repeat{
    coverage <- c(coverage, mean(pileupSubset[c(X:Y), 2]))
    X <- X + (windowSize / 100)
    Y <- Y + (windowSize / 100)
    if (Y > nrow(pileupSubset)) break
  }
  if (mode == "genome" || chunkContigs == TRUE) {
    position <- seq(pileupSubset[1, 3], pileupSubset[nrow(pileupSubset), 3], length.out = length(coverage))
  } else {
    position <- seq(windowSize, length(coverage) * windowSize, windowSize)
  }
  refName <- rep(pileupSubset[1, 1], length(position))
  newdataset <- cbind.data.frame(refName, coverage, position)
  newdataset[do.call(cbind, lapply(newdataset, is.nan))] <- 0
  return(newdataset)
}

#' Reformat input pileup file
#'
#' Place columns in correct order, clean accessions by removing text after white
#' space, and name columns
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
pileupFormatter <- function(pileup, mode) {
  colClasses <-
    lapply(seq_along(pileup), function(i) {
      class(pileup[, i])
    })
  for (i in c(which(colClasses == "integer"))) {
    if (mode == "genome") {
      if (pileup[1, i] == 0) {
        posColIdx <- i
      }
    } else {
      if (length(which(pileup[, i] == 100)) > 1) {
        posColIdx <- i
      }
    }
  }
  cleanPileup <-
    cbind.data.frame(
      pileup[, which(colClasses == "character")],
      pileup[, which(colClasses == "numeric")],
      pileup[, posColIdx]
    )
  colnames(cleanPileup) <- c("refName", "coverage", "position")
  cleanPileup$refName <- gsub("\\s.*", "", cleanPileup$refName)
  return(cleanPileup)
}

