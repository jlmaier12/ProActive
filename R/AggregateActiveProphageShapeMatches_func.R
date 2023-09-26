#' Collects prophage prediction shape-match information
#'
#' Collects shape information associated with all contigs predicted as containing prophages based on the results from the shape_matcher function.
#'
#' @param best_match_list Predictions made with shape_matcher function. Predictions are stored as the first item in the best_match_list.
#'
#' @keywords internal
allprophages_func_WC <- function(best_match_list){
  A<-1
  prophageprediction_list <- list()
  for (i in seq(1,length(best_match_list),1)){
    prediction <-  best_match_list[[i]][[6]]
    if(prediction=="None") next
    prophageprediction_list[[A]] <- best_match_list[[i]]
    A <- A+1
  }
  return(prophageprediction_list)
}
