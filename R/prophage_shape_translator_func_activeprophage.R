#' Full prophage shape shape-translator for whole-community read coverage only
#'
#' Translates full prophage shape across a contig. Translate the shape 2000bp at a time. Stop translating when the prophage shape is 5000bp from the end of the contig.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param best_match_info The information associated with the current best shape match. Includes the comparison score, the minimum and maximum shape y-acis values, the shape length, and the shape start and stop positions on the contig
#' @param windowsize The window size used to re-average read coverage datasets
#' @param shape A vector containing the values associated with the prophage shape to be translated across the contig
proshape_translator_func_WC <- function(microbial_subset, best_match_info, windowsize, shape){
  Cov_values_contig <- microbial_subset[,2]
  min_shape_cov <- min(shape)
  max_shape_cov <- max(shape)
  repeat {
    shape <- c(rep(min_shape_cov, (2000/windowsize)),shape[-c((length(shape)-((2000/windowsize)-1)):length(shape))])
    if(shape[length(shape)-(15000/windowsize)]>min_shape_cov) break
    diff <- mean(abs(Cov_values_contig - shape))
    start_pos <- which(shape==max(shape))[1]
    end_pos <- which(shape==max(shape))[length(which(shape==max(shape)))]
    if (diff < best_match_info[[1]]){
      best_match_info <- list(diff, min_shape_cov, max_shape_cov, start_pos, end_pos)
    }
  }
  return(best_match_info)
}
