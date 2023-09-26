#' Prophage off-left shape-matcher for whole-community read coverage
#'
#' Build and translate and prophage shape going off the left side of the graph across a contig. Stop translating when the shape left on the contig is 10,000bp. Translate the shape 2000bp at a time.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#'
#' @keywords internal
prophage_off_left_func_WC <- function (microbial_subset, windowsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  threequarter_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  quarter_read_cov <- (abs(max_read_cov - min_read_cov))/4
  maxread_cov_steps <- (abs(max_read_cov - (min_read_cov+threequarter_read_cov)))/5
  minread_cov_steps <- (abs((min_read_cov+quarter_read_cov) - min_read_cov))/5
  Cov_values_contig <- microbial_subset[,2]
  startingcoverages <- seq((min_read_cov+threequarter_read_cov), max_read_cov, maxread_cov_steps)
  startingmincoverages <- seq(min_read_cov,(min_read_cov+quarter_read_cov), minread_cov_steps)
  shape <- c(rep(startingcoverages[1], nrow(microbial_subset)-(10000/windowsize)), rep(min_read_cov, 10000/windowsize))
  diff <- mean(abs(Cov_values_contig - shape))
  end_pos <- (which(shape == min(shape))[1])-1
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], 1, end_pos)
  for(cov in startingcoverages) {
    min_read_cov <- min(microbial_subset[,2])
    shape <- c(rep(cov, nrow(microbial_subset)-(10000/windowsize)), rep(min_read_cov, 10000/windowsize))
    repeat {
      if (length(which(shape==cov)) < 10000/windowsize) break
      if (diff < best_match_info[[1]]){
        best_match_info <- list(diff, min_read_cov, cov, 1, end_pos)
      }
      shape <- c(shape[-c(1:(2000/windowsize))], rep(min_read_cov, (2000/windowsize)))
      diff <- mean(abs(Cov_values_contig - shape))
      end_pos <- (which(shape == min_read_cov)[1])-1
    }
    for(mincov in startingmincoverages) {
      shape <- c(rep(cov, nrow(microbial_subset)-(10000/windowsize)), rep(mincov, 10000/windowsize))
      repeat {
        if (length(which(shape==cov)) < 10000/windowsize) break
        if (diff < best_match_info[[1]]){
          best_match_info <- list(diff, mincov, cov, 1, end_pos)
        }
        shape <- c(shape[-c(1:(2000/windowsize))], rep(mincov, (2000/windowsize)))
        diff <- mean(abs(Cov_values_contig - shape))
        end_pos <- (which(shape == mincov)[1])-1
      }
    }
  }
  best_match_results <- append(best_match_info, "prophage")
  return(best_match_results)
}
