#' Change the pileup window size
#'
#' Re-averages windows of pileup files with 100bp windows to reduce pileup size.
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
changewindowSize <- function(pileup, windowSize, mode) {
  coverage <- vector()
  X <- 0
  Y <- windowSize / 100
  repeat{
    coverage <- c(coverage, mean(pileup[c(X:Y), 2]))
    X <- X + (windowSize / 100)
    Y <- Y + (windowSize / 100)
    if (Y > nrow(pileup)) break
  }
  if (mode == "genome") {
    position <- seq(pileup[1, 3], pileup[nrow(pileup), 3], length.out = length(coverage))
  } else {
    position <- seq(windowSize, length(coverage) * windowSize, windowSize)
  }
  refName <- rep(pileup[1, 1], length(position))
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
        posColIdx <<- i
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

