#' shape-match size calculator
#'
#' Calculate the size, in base pairs, of the matching region for shapes detected in whole-community read coverages
#'
#' @param predictionsummary Prediction summary table
#' @param prophageprediction_list A list containing shape information associated with all contigs containing a potential active prophages. Generated with the allprophages_func
#' @param windowsize The window size used to re-average read coverage datasets
#'
#' @keywords internal
matchsize_checker <- function(predictionsummary, prophageprediction_list, windowsize){
  ref_name <- rep(NA, length(prophageprediction_list))
  match_size <- rep(NA, length(prophageprediction_list))
  for (i in seq(1,length(prophageprediction_list),1)) {
    ref_name[i] <- prophageprediction_list[[i]][[8]]
    start_pos <- prophageprediction_list[[i]][[4]]
    end_pos <- prophageprediction_list[[i]][[5]]
    match_size[i] <- (length(c(start_pos:end_pos))-1) *windowsize
  }
  match_lengthsummary_table <- cbind.data.frame(match_size, ref_name)
  summary_allcontigs_withmatchlength <- merge(predictionsummary, match_lengthsummary_table, by="ref_name", all.x=TRUE)
  return(summary_allcontigs_withmatchlength)
}


