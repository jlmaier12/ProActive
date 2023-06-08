#' No transduction shape-matcher for whole-community read coverage
#'
#' Assess whether a contig does not have a read coverage pattern associated with a transduction event. A horizontal line at the mean coverage should be an optimal match if the contig read coverage displays no transduction patterns
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
notransduction_func_WC <- function (microbial_subset) {
  shape <- rep(mean(microbial_subset[,2]),nrow(microbial_subset))
  Cov_values_contig <- microbial_subset[,2]
  best_match_score <- mean(abs(Cov_values_contig - shape))
  best_cov <- mean(microbial_subset[,2])
  for (cov in seq(min(microbial_subset[,2]), max(microbial_subset[,2]),4))
    shape <- rep(cov,nrow(microbial_subset))
  Cov_values_contig <- microbial_subset[,2]
  diff <- mean(abs(Cov_values_contig - shape))
  if (diff < best_match_score) {
    best_match_score <- diff
    best_cov <- cov
  }
  best_match_info <- list(best_match_score, best_cov, nrow(microbial_subset), "NA", "NA", "None")
  return(best_match_info)
}
