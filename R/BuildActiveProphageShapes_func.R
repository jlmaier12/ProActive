#' Shape-builder for active prophages detected in whole-community read coverages
#'
#' Builds the shape (vector) associated with the 'best shape-match' information associated with each contig so it can be plotted.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param prophagepredictions A list containing shape information associated with all contigs containing a potential active prophage. Generated with the allprophages_func
#' @param i The index position associated with the current contig's best shape-match information
shape_builder_func_WC <- function(microbial_subset, prophagepredictions, i){
  min_read_cov <- prophagepredictions[[i]][[2]]
  max_read_cov <- prophagepredictions[[i]][[3]]
  start_pos <- prophagepredictions[[i]][[4]]
  end_pos <- prophagepredictions[[i]][[5]]
  prediction <- prophagepredictions[[i]][[6]]
  if (start_pos==1) {
    shape <- c(rep(max_read_cov,end_pos), rep(min_read_cov, (nrow(microbial_subset)-end_pos)))
  } else if (end_pos == nrow(microbial_subset)){
    shape <- c(rep(min_read_cov, start_pos), rep(max_read_cov, (nrow(microbial_subset)-start_pos)))
  } else{
    match_region <- end_pos-start_pos
    shape <- c(rep(min_read_cov, start_pos), rep(max_read_cov, match_region), rep(min_read_cov, (nrow(microbial_subset)-(match_region+start_pos))))
  }
  return(shape)
}
