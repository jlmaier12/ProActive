#' Summarizes best shape-matches into a table
#'
#' Summarizes the predictions made in the shape_matcher function based on which shape match achieved the lowest mean absolute difference for each contig. Outputs results in a table rather than a list.
#'
#' @param best_match_list Predictions made with shape_matcher function. Predictions are stored as the first item in the best_match_list.
#' @param metagenome_pileup A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The number of base pairs to average coverage values over. For compatibility with shape-matching, window sizes can only be 200, 500, 1000 and 2000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowsize is the default.
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
contig_prediction_summary_WC <- function(metagenome_pileup,best_match_list, windowsize, mode){
  ref_name <- rep(NA, length(best_match_list))
  confidence <- rep(NA, length(best_match_list))
  elev_ratio <- rep(NA, length(best_match_list))
  start_pos <- rep(NA, length(best_match_list))
  stop_pos <- rep(NA, length(best_match_list))
  prediction <- rep(NA, length(best_match_list))
  if(length(best_match_list)==0){
    cat("no elevated read coverage detected \n")
  }
  for (i in seq(1,length(best_match_list),1)){
    ref_name[i] <- best_match_list[[i]][[8]]
    contig_name <- best_match_list[[i]][[8]]
    contig <- metagenome_pileup[which(metagenome_pileup[,1]==contig_name),]
    contig <- windowsize_func(contig, windowsize, mode)
    confidence[i] <- best_match_list[[i]][[9]]
    elev_ratio[i] <- best_match_list[[i]][[6]]
    prediction[i] <- best_match_list[[i]][[7]]
    start_pos[i] <- contig[best_match_list[[i]][[4]],3]
    stop_pos[i] <- contig[best_match_list[[i]][[5]],3]
  }
  Prediction_summary <- cbind.data.frame(ref_name, prediction, confidence, elev_ratio, start_pos, stop_pos)
  return(Prediction_summary)
}
