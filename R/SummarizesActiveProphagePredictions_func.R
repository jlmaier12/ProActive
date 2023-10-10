#' Summarizes best shape-matches into a table
#'
#' Summarizes the predictions made in the shape_matcher function based on which shape match achieved the lowest mean absolute difference for each contig. Outputs results in a table rather than a list.
#'
#' @param best_match_list Predictions made with shape_matcher function. Predictions are stored as the first item in the best_match_list.
#'
#' @keywords internal
contig_prediction_summary_WC <- function(best_match_list){
  ref_name <- rep(NA, length(best_match_list))
  confidence <- rep(NA, length(best_match_list))
  elev_ratio <- rep(NA, length(best_match_list))
  if(length(best_match_list)==0){
    print("no elevated read coverage detected")
  }
  for (i in seq(1,length(best_match_list),1)){
    ref_name[i] <- best_match_list[[i]][[8]]
    confidence[i] <- best_match_list[[i]][[9]]
    elev_ratio[i] <- best_match_list[[i]][[6]]
  }
  Prediction_summary <- cbind(ref_name, confidence, elev_ratio)
  return(Prediction_summary)
}
