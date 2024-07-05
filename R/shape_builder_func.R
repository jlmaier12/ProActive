#' Shape-builder for patterns detected in whole-community read coverages
#'
#' Builds the shape (vector) associated with the 'best shape-match' information associated with each contig so it can be plotted.
#'
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param ProActivePredictions A list containing shape information associated with all contigs classified by ProActive
#' @param i The index position associated with the current contig's best shape-match information
#'
#' @keywords internal
shape_builder_func <- function(microbial_subset, ProActivePredictions, i){
  min_read_cov <- ProActivePredictions[[i]][[2]]
  max_read_cov <- ProActivePredictions[[i]][[3]]
  start_pos <- ProActivePredictions[[i]][[4]]
  end_pos <- ProActivePredictions[[i]][[5]]
  prediction <- ProActivePredictions[[i]][[7]]
  if(prediction == "Gap") {
    if (start_pos==1) {
      shape <- c(rep(min_read_cov,end_pos), rep(max_read_cov, (nrow(microbial_subset)-end_pos)))
    } else if (end_pos == nrow(microbial_subset)){
      shape <- c(rep(max_read_cov, start_pos), rep(min_read_cov, (nrow(microbial_subset)-start_pos)))
    } else{
      match_region <- end_pos-start_pos
      shape <- c(rep(max_read_cov, start_pos), rep(min_read_cov, match_region), rep(max_read_cov, (nrow(microbial_subset)-(match_region+start_pos))))
    }
  } else if (prediction == "None"){
      shape <- rep(min_read_cov, nrow(microbial_subset))
  } else {
  if (start_pos==1) {
    shape <- c(rep(max_read_cov,end_pos), rep(min_read_cov, (nrow(microbial_subset)-end_pos)))
  } else if (end_pos == nrow(microbial_subset)){
    shape <- c(rep(min_read_cov, start_pos), rep(max_read_cov, (nrow(microbial_subset)-start_pos)))
  } else{
    match_region <- end_pos-start_pos
    shape <- c(rep(min_read_cov, start_pos), rep(max_read_cov, match_region), rep(min_read_cov, (nrow(microbial_subset)-(match_region+start_pos))))
  }
  }
  return(shape)
}
