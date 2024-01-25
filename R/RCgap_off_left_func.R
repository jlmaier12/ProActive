#' Read coverage gap off-left shape-matcher for whole-community read coverage
#'
#' Build and translate and read coverage gap shape going off the left side of the graph across a contig. Stop translating when the shape left on the contig is 10,000bp. Translate the shape 2000bp at a time.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minsize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxsize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#'
#' @keywords internal
RCgap_off_left_func <- function (microbial_subset, windowsize, minsize, maxsize) {
  max_read_cov <- max(microbial_subset[,2])
  min_read_cov <- min(microbial_subset[,2])
  threequarter_read_cov <- (abs(max_read_cov-min_read_cov)*3)/4
  quarter_read_cov <- (abs(max_read_cov - min_read_cov))/4
  maxread_cov_steps <- (abs(max_read_cov - (min_read_cov+threequarter_read_cov)))/5
  minread_cov_steps <- (abs((min_read_cov+quarter_read_cov) - min_read_cov))/5
  Cov_values_contig <- microbial_subset[,2]
  startingcoverages <- seq((min_read_cov+threequarter_read_cov), max_read_cov, maxread_cov_steps)
  startingmincoverages <- seq(min_read_cov,(min_read_cov+quarter_read_cov), minread_cov_steps)
  shape_length <- ifelse ((nrow(microbial_subset)-(10000/windowsize))>(maxsize/windowsize),maxsize/windowsize,nrow(microbial_subset)-(10000/windowsize))
  nonshape <- nrow(microbial_subset)-shape_length
  shape <- c(rep(min_read_cov, shape_length), rep(startingcoverages[1], nonshape))
  diff <- mean(abs(Cov_values_contig - shape))
  end_pos <- (which(shape == max(shape))[1])-1
  elevation_ratio <- max(shape)/min(shape)
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], 1, end_pos, elevation_ratio)
  for(cov in startingcoverages) {
    min_read_cov <- min(microbial_subset[,2])
    shape <- c(rep(min_read_cov, shape_length),rep(cov, nonshape))
    repeat {
      if (diff < best_match_info[[1]]){
        elevation_ratio <- max(shape)/min(shape)
        best_match_info <- list(diff, min_read_cov, cov, 1, end_pos, elevation_ratio)
      }
      if (length(which(shape==min_read_cov))<= 2000/windowsize) break
      shape <- c(shape[-c(1:(2000/windowsize))], rep(cov, (2000/windowsize)))
      if (length(which(shape==min_read_cov)) < minsize/windowsize) break
      diff <- mean(abs(Cov_values_contig - shape))
      end_pos <- (which(shape == cov)[1])-1
    }
    for(mincov in startingmincoverages){
      shape <- c(rep(mincov, shape_length), rep(cov, nonshape))
      repeat {
        if (diff < best_match_info[[1]]){
          elevation_ratio <- max(shape)/min(shape)
          best_match_info <- list(diff, mincov, cov, 1, end_pos, elevation_ratio)
        }
        if (length(which(shape==mincov))<= 2000/windowsize) break
        shape <- c(shape[-c(1:(2000/windowsize))], rep(cov, (2000/windowsize)))
        if (length(which(shape==mincov)) < minsize/windowsize) break
        diff <- mean(abs(Cov_values_contig - shape))
        end_pos <- (which(shape == cov)[1])-1
      }
    }
  }
  best_match_results <- c(best_match_info, "Gap")
  return(best_match_results)
}
