#' Determines pattern-match quality
#'
#' The quality is either 'VeryHigh', 'High', or 'Low'. Pattern-matches
#' receive a very high quality rating when all the lowest or highest coverage values on
#' the contig/chunk fall with the gapped or elevated region, respectively, of the
#' pattern-match. Pattern-matches receive a high quality rating when the contig/chunk's
#' maximum or minimum read coverage value falls within the elevated or gapped region,
#' respectively, of the pattern-match. Pattern-matches receive a low quality rating
#' when the contig/chunk's maximum or minimum read coverage value does not fall within
#' the elevated or gapped region, respectively, of the pattern-match.
#'
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param mode Either "genome" or "metagenome".
#' @keywords internal
classifQuality <- function(bestMatchList, pileup, windowSize, mode) {
  lapply(seq_along(bestMatchList), function(i) {
    bestMatchInfo <- bestMatchList[[i]]
    classification <- bestMatchInfo[[7]]
    refName <- bestMatchInfo[[8]]
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    matchRegion <- c(bestMatchInfo[[4]]:bestMatchInfo[[5]])
    if (classification == "Gap") {
      covIndex <- which(pileupSubset[, 2] == min(pileupSubset[, 2]))[1]
      covs <- order(pileupSubset[, 2], decreasing = FALSE)[1:length(matchRegion)]
    } else if (classification == "Elevation") {
      covIndex <- which(pileupSubset[, 2] == max(pileupSubset[, 2]))[1]
      covs <- order(pileupSubset[, 2], decreasing = TRUE)[1:length(matchRegion)]
    } else {
      bestMatchList[[i]] <<- c(bestMatchInfo, NA)
      return(NULL)
    }
    if (FALSE %in% (covs %in% matchRegion) == FALSE) {
      bestMatchList[[i]] <<- c(bestMatchInfo, "VeryHigh")
    } else if ((covIndex %in% matchRegion) == TRUE) {
      bestMatchList[[i]] <<- c(bestMatchInfo, "High")
    } else {
      bestMatchList[[i]] <<- c(bestMatchInfo, "Low")
    }
  })
  return(bestMatchList)
}

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
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
classifSumm <- function(pileup, bestMatchList, windowSize, mode) {
  refName <- rep(NA, length(bestMatchList))
  confidence <- rep(NA, length(bestMatchList))
  elevRatio <- rep(NA, length(bestMatchList))
  startPos <- rep(NA, length(bestMatchList))
  endPos <- rep(NA, length(bestMatchList))
  classification <- rep(NA, length(bestMatchList))
  matchSize <- rep(NA, length(bestMatchList))
  if (length(bestMatchList) == 0) {
    stop("No pattern-matches detected")
  }
  lapply(seq_along(bestMatchList), function(i) {
    refName[i] <<- bestMatchList[[i]][[8]]
    pileupSubset <- pileup[which(pileup[, 1] == bestMatchList[[i]][[8]]), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    confidence[i] <<- bestMatchList[[i]][[9]]
    elevRatio[i] <<- bestMatchList[[i]][[6]]
    classification[i] <<- bestMatchList[[i]][[7]]
    startPos[i] <<- pileupSubset[bestMatchList[[i]][[4]], 3]
    endPos[i] <<- pileupSubset[bestMatchList[[i]][[5]], 3]
    matchSize[i] <<- (length(seq(pileupSubset[bestMatchList[[i]][[4]], 3], pileupSubset[bestMatchList[[i]][[5]], 3], windowSize)) - 1) * windowSize
  })
  classifSummTable <- cbind.data.frame(refName, classification, confidence, elevRatio, startPos, endPos, matchSize)
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
  newBestMatchList <- list()
  length(bestMatchList)
  X <- 1
  lapply(seq_along(bestMatchList), function(i) {
    bestMatchInfo <- bestMatchList[[i]]
    classification <- bestMatchInfo[[7]]
    if (classification == "NoPattern") {
      return(NULL)
    } else {
      newBestMatchList[[X]] <<- bestMatchInfo
      X <<- X + 1
    }
  })

  return(newBestMatchList)
}
