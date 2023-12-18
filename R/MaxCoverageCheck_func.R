#' Checks for max coverage outside match region
#'
#' @param best_match_list
#' @param microbial_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The number of base pairs to average coverage values over. For compatibility with shape-matching, window sizes can only be 200, 500, 1000 and 2000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowsize is the default.
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
max_coverage_check <- function (best_match_list, microbialread_dataset, windowsize, mode){
  ProActiveConfPredictions <- list()
  ProActivePredictions <- list()
  ProActiveHonMentions <- list()
  NonePredictions <- list()
  A<-1
  B<-1
  C<-1
  D<-1
  for(index in seq(1,length(best_match_list),1)) {
    best_match_info <- best_match_list[[index]]
    prediction <- best_match_info[[7]]
    if(prediction=="None") {
    NonePredictions[[C]] <- best_match_info
    C <- C+1
    next
    }
    ref_name <- best_match_info[[8]]
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1]==ref_name),]
    microbial_subset <- windowsize_func(microbial_subset, windowsize, mode)
    max_coverage_index <- which(microbial_subset[,2]==max(microbial_subset[,2]))[1]
    match_region <- c(best_match_info[[4]]:best_match_info[[5]])
    max_covs <- order(microbial_subset[,2], decreasing=TRUE)[1:length(match_region)]
    elevation_ratio <- best_match_info[[6]]
    if (FALSE %in% (max_covs %in% match_region) == FALSE){
      ProActiveConfPredictions[[D]] <- c(best_match_info, "VeryHigh")
      D <- D+1
    }else if((max_coverage_index %in% match_region) == TRUE){
      ProActivePredictions[[A]] <- c(best_match_info,"High")
      A <- A+1
    }else{
      ProActiveHonMentions[[B]] <- c(best_match_info, "Low")
      B <- B+1
    }
  }
    max_cov_check_results <- list(ProActiveConfPredictions, ProActivePredictions,ProActiveHonMentions,NonePredictions)
    return(max_cov_check_results)
}



