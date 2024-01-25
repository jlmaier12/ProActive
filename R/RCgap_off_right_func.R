#' Read coverage gap off-right shape-matcher
#'
#' Build and translate read coverage gap shape going off the right side of the graph across a contig. Stop translating when the shape left on the contig is 10,000bp. Translate the shape 2000bp at a time.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minsize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxsize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#'
#' @keywords internal
RCgap_off_right_func <- function (microbial_subset, windowsize, minsize, maxsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  threequarter_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  quarter_read_cov <- (abs(max_read_cov - min_read_cov))/4
  maxread_cov_steps <- (abs(max_read_cov - (min_read_cov+threequarter_read_cov)))/5
  minread_cov_steps <- (abs((min_read_cov+quarter_read_cov) - min_read_cov))/5
  Cov_values_contig <- microbial_subset[,2]
  startingcoverages <- seq((min_read_cov+threequarter_read_cov), max_read_cov, maxread_cov_steps)
  startingmincoverages <- seq(min_read_cov,(min_read_cov+quarter_read_cov), minread_cov_steps)
  shape_length <- ifelse((nrow(microbial_subset)-(10000/windowsize))>maxsize/windowsize, maxsize/windowsize, nrow(microbial_subset)-(10000/windowsize))
  nonshape <- nrow(microbial_subset)-shape_length
  shape <- c(rep(startingcoverages[1], nonshape),rep(min_read_cov, shape_length))
  diff <- mean(abs(Cov_values_contig - shape))
  start_pos <- (which(shape == min(shape))[1])
  elevation_ratio <- max(shape)/min(shape)
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], start_pos, length(shape), elevation_ratio)
  for(cov in startingcoverages) {
    min_read_cov <- min(microbial_subset[,2])
    shape <- c(rep(cov, nonshape), rep(min_read_cov, shape_length))
    repeat {
      if (diff < best_match_info[[1]]){
        elevation_ratio <- max(shape)/min(shape)
        best_match_info <- list(diff, min_read_cov, cov, start_pos, length(shape), elevation_ratio)
      }
      if (length(which(shape==min_read_cov))<= 2000/windowsize) break
      shape <- c(rep(cov,(2000/windowsize)),shape[-c(((length(shape))-((2000/windowsize)-1)):length(shape))]) #variable, removing 2000bp at a time
      if (length(which(shape==min_read_cov)) < minsize/windowsize) break
      diff <- mean(abs(Cov_values_contig - shape))
      start_pos <- (which(shape == min(shape))[1])
    }
    for(mincov in startingmincoverages) {
      shape <- c(rep(cov, nonshape), rep(mincov, shape_length))
      repeat {
        if (diff < best_match_info[[1]]){
          elevation_ratio <- max(shape)/min(shape)
          best_match_info <- list(diff, mincov, cov, start_pos, length(shape), elevation_ratio)
        }
        if (length(which(shape==mincov))<= 2000/windowsize) break
        shape <- c(rep(cov,(2000/windowsize)),shape[-c(((length(shape))-((2000/windowsize)-1)):length(shape))]) #variable, removing 2000bp at a time
        if (length(which(shape==mincov)) < minsize/windowsize) break
        diff <- mean(abs(Cov_values_contig - shape))
        start_pos <- (which(shape == min(shape))[1])
      }
    }
  }
  best_match_results <- c(best_match_info, "Gap")
  return(best_match_results)
}
