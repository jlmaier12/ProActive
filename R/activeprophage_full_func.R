#' Full prophage shape-matcher
#'
#' Build the full prophage shape and provide it as input to the proshape_translator_func to be translated across a contig. The shape is made smaller length-wise 2000bp at a time, and the shape stops decreasing in length once it reaches 10000bp.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
full_prophage_func_WC <- function (microbial_subset, windowsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  half_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  bottomtotop_read_cov <- abs(max_read_cov - half_read_cov)/4
  Cov_values_contig <- microbial_subset[,2]
  startingcoverages <- seq(half_read_cov, max_read_cov, bottomtotop_read_cov)
  shape <- c(rep(min_read_cov, 15000/windowsize), rep(startingcoverages[1], nrow(microbial_subset)-(15000/windowsize)))
  diff <- mean(abs(Cov_values_contig - shape))
  start_pos <- (which(shape == max(shape))[1])
  end_pos <- which(shape==max(shape))[length(which(shape==max(shape)))]
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], start_pos, end_pos)
  for(cov in startingcoverages) {
    shape <- c(rep(min_read_cov, 15000/windowsize), rep(cov, nrow(microbial_subset)-(15000/windowsize)))
    repeat {
      best_match_info <- proshape_translator_func_WC(microbial_subset, best_match_info, windowsize, shape)
      middle_rows <- which(shape == cov)
      shape <- c(shape[-c(middle_rows[2]:middle_rows[(2000/windowsize)+1])],rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
      if (length(middle_rows) < 5000/windowsize) break
    }
  }
  best_match_results <- append(best_match_info, "Prophage")
  return(best_match_results)
}
