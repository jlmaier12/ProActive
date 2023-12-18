#' Full prophage shape-matcher
#'
#' Build the full prophage shape and provide it as input to the proshape_translator_func to be translated across a contig. The shape is made smaller length-wise 2000bp at a time, and the shape stops decreasing in length once it reaches 10000bp.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minsize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxsize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#'
#' @keywords internal
full_prophage_func_WC <- function (microbial_subset, windowsize, minsize, maxsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  threequarter_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  quarter_read_cov <- (abs(max_read_cov - min_read_cov))/4
  maxread_cov_steps <- (abs(max_read_cov - (min_read_cov+threequarter_read_cov)))/5
  minread_cov_steps <- (abs((min_read_cov+quarter_read_cov) - min_read_cov))/5
  Cov_values_contig <- microbial_subset[,2]
  startingcoverages <- seq((min_read_cov+threequarter_read_cov), max_read_cov, maxread_cov_steps)
  startingmincoverages <- seq(min_read_cov,(min_read_cov+quarter_read_cov), minread_cov_steps)
  shape_length <- ifelse((nrow(microbial_subset)-(20000/windowsize))>(maxsize/windowsize), maxsize/windowsize, nrow(microbial_subset)-(20000/windowsize))
  nonshape <- nrow(microbial_subset)-(shape_length+(10000/windowsize))
  shape <- c(rep(min_read_cov, 10000/windowsize), rep(startingcoverages[1], shape_length), rep(min_read_cov, nonshape))
  diff <- mean(abs(Cov_values_contig - shape))
  start_pos <- (which(shape == max(shape))[1])
  end_pos <- which(shape==max(shape))[length(which(shape==max(shape)))]
  elevation_ratio <- max(shape)/min(shape)
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], start_pos, end_pos, elevation_ratio)
  for(cov in startingcoverages) {
    min_read_cov <- min(microbial_subset[,2])
    shape <- c(rep(min_read_cov, 10000/windowsize), rep(cov, shape_length), rep(min_read_cov, nonshape))
    repeat {
      middle_rows <- which(shape == cov)
      if (length(middle_rows) < minsize/windowsize) break
      best_match_info <- proshape_translator_func_WC(Cov_values_contig, best_match_info, windowsize, shape)
      if (length(middle_rows)<= 2000/windowsize) break
      shape <- c(shape[-c(middle_rows[2]:middle_rows[(2000/windowsize)+1])], rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
    }
    for(mincov in startingmincoverages) {
      shape <- c(rep(mincov, 10000/windowsize), rep(cov, shape_length), rep(mincov, nonshape))
      repeat {
        middle_rows <- which(shape == cov)
        if (length(middle_rows) < minsize/windowsize) break
        best_match_info <- proshape_translator_func_WC(Cov_values_contig, best_match_info, windowsize, shape)
        if (length(middle_rows)<= 2000/windowsize) break
        shape <- c(shape[-c(middle_rows[2]:middle_rows[(2000/windowsize)+1])], rep(mincov,2000/windowsize)) #remove 2000bp at a time
      }
    }
  }
  best_match_results <- c(best_match_info, "prophage")
  return(best_match_results)
}
