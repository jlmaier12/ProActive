#' No read coverage pattern
#'
#' Assess whether a contig/chunk does not have an elevated/gapped read coverage pattern.
#' A horizontal line at the mean or median coverage should be the best
#' match if the contig/chunk read coverage is not gapped or elevated.
#'
#' @param pileupSubset A subset of the read coverage dataset that pertains only to
#' the contig currently being assessed
#' @importFrom stats median
#' @keywords internal
noPattern <- function(pileupSubset) {
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
