#' No transduction pattern-matcher for whole-community read coverage
#'
#' Assess whether a contig does not have a read coverage pattern associated with a transduction event. A horizontal line at the mean coverage should be an optimal match if the contig read coverage displays no transduction patterns
#'
#' @param pileupSubset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#'
#' @keywords internal
noPattern <- function (pileupSubset) {
    pattern1 <- rep(median(pileupSubset[, 2]), nrow(pileupSubset))
    pattern2 <- rep(mean(pileupSubset[, 2]), nrow(pileupSubset))
    diff1 <- mean(abs(pileupSubset[, 2] - pattern1))
    diff2 <- mean(abs(pileupSubset[, 2] - pattern2))
    diff <- ifelse(diff1 < diff2, diff1, diff2)
    value <- ifelse(diff1 < diff2, median(pileupSubset[, 2]), mean(pileupSubset[, 2]))
    bestMatchInfo <-
      list(
        diff,
        value,
        nrow(pileupSubset),
        1,
        nrow(pileupSubset),
        NA,
        "NoPattern"
      )
    return(bestMatchInfo)
  }


