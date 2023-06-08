#' Prophage off-right shape-matcher
#'
#' Build and translate and prophage shape going off the right side of the graph across a contig. Stop translating when the shape left on the contig is 10,000bp. Translate the shape 2000bp at a time.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
prophage_off_right_func_WC <- function (microbial_subset, windowsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  half_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  bottomtotop_read_cov <- abs(max_read_cov - half_read_cov)/4
  Cov_values_contig <- microbial_subset[,2]
  if (abs(half_read_cov-min_read_cov)<=10) {
    half_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  }
  startingcoverages <- seq(half_read_cov,  max_read_cov, bottomtotop_read_cov)
  shape <- c(rep(min_read_cov, 15000/windowsize), rep(startingcoverages[1], nrow(microbial_subset)-(15000/windowsize)))
  diff <- mean(abs(Cov_values_contig - shape))
  start_pos <- (which(shape == max(shape))[1])
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], start_pos, length(shape))
  for(cov in startingcoverages) {
    shape <- c(rep(min_read_cov, 15000/windowsize), rep(cov, nrow(microbial_subset)-(15000/windowsize)))
    repeat {
      if (length(which(shape==cov)) < 10000/windowsize) break
      if (diff < best_match_info[[1]]){
        best_match_info <- list(diff, min_read_cov, cov, start_pos, length(shape))
      }
      shape <- c(rep(min_read_cov,(2000/windowsize)),shape[-c(((length(shape))-((2000/windowsize)-1)):length(shape))]) #variable, removing 2000bp at a time
      diff <- mean(abs(Cov_values_contig - shape))
      start_pos <- (which(shape == max(shape))[1])
    }
  }
  best_match_results <- append(best_match_info, "Prophage")
  return(best_match_results)
}
