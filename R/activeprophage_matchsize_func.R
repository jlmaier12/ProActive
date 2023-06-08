#' Active prophage shape-match size calculator
#'
#' Calculate the size, in base pairs, of the matching region for prophage shapes detected in whole-community read coverages
#'
#' @param predictionsummary Prediction summary table
#' @param prophageprediction_list A list containing shape information associated with all contigs containing a potential active prophages. Generated with the allprophages_func
#' @param windowsize The window size used to re-average read coverage datasets
activeprophage_matchsize_checker <- function(predictionsummary, prophageprediction_list, windowsize){
  ref_name <- rep(NA, length(prophageprediction_list))
  match_size <- rep(NA, length(prophageprediction_list))
  for (i in seq(1,length(prophageprediction_list),1)) {
    ref_name[i] <- prophageprediction_list[[i]][[7]]
    start_pos <- prophageprediction_list[[i]][[4]]
    end_pos <- prophageprediction_list[[i]][[5]]
    match_size[i] <- (end_pos-start_pos) *windowsize
  }
  match_lengthsummary_table <- cbind(match_size, ref_name)
  summary_allcontigs_withmatchlength <- merge(predictionsummary, match_lengthsummary_table, by="ref_name", all.x=TRUE)
  return(summary_allcontigs_withmatchlength)
}
