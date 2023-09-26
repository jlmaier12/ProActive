#' Checks for max coverage outside match region
#'
#' @param best_match_list
#' @param microbial_subset
#'
#' @keywords internal
max_coverage_check <- function (best_match_list, microbialread_dataset, windowsize){
  ProActivePredictions <- list()
  ProActiveHonMentions <- list()
  NonePredictions <- list()
  A<-1
  B<-1
  C<-1
  for(index in seq(1,length(best_match_list),1)) {
    best_match_info <- best_match_list[[index]]
    prediction <- best_match_info[[6]]
    if(prediction=="None") {
    NonePredictions[[C]] <- best_match_info
    C <- C+1
    next
    }
    ref_name <- best_match_info[[7]]
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1]==ref_name),]
    microbial_subset <- windowsize_func(microbial_subset, windowsize)
    max_coverage_index <- which(microbial_subset[,2]==max(microbial_subset[,2]))[1]
    match_region <- c(best_match_info[[4]]:best_match_info[[5]])
    if((max_coverage_index %in% match_region) == TRUE){
      ProActivePredictions[[A]] <- append(best_match_info,"Yes")
      A <- A+1
    } else{
      ProActiveHonMentions[[B]] <- append(best_match_info, "No")
      B <- B+1
    }
  }
    max_cov_check_results <- list(ProActivePredictions,ProActiveHonMentions,NonePredictions)
    return(max_cov_check_results)
}




