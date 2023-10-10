#' No transduction shape-matcher for whole-community read coverage
#'
#' Assess whether a contig does not have a read coverage pattern associated with a transduction event. A horizontal line at the mean coverage should be an optimal match if the contig read coverage displays no transduction patterns
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#'
#' @keywords internal
notransduction_func_WC <- function (microbial_subset) {
  shape <- rep(mean(microbial_subset[,2]),nrow(microbial_subset))
  Cov_values_contig <- microbial_subset[,2]
  best_match_score <- mean(abs(Cov_values_contig - shape))
  best_cov <- mean(microbial_subset[,2])
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  quarter_read_cov <- (abs(max_read_cov - min_read_cov))/4
  minread_cov_steps <- (abs((min_read_cov+quarter_read_cov) - min_read_cov))/5
  startingmincoverages <- seq(min_read_cov,(min_read_cov+quarter_read_cov), minread_cov_steps)
  for (cov in startingmincoverages)
    shape <- rep(cov,nrow(microbial_subset))
    diff <- mean(abs(Cov_values_contig - shape))
  if (diff < best_match_score) {
    best_match_score <- diff
    best_cov <- cov
  }
  best_match_info <- list(best_match_score, best_cov, nrow(microbial_subset), "NA", "NA", "NA", "None")
  return(best_match_info)
}
