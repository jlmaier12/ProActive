#' Shape-matching function for active MGE in mapped read coverages
#'
#' Creates the microbial_subset, representative of one contig, that is used as input for each individual shape-building function. After the information associated with the best match for each shape is obtained, the shape with the lowest mean absolute difference (score) is chosen as the prediction for the contig being assessed.
#'
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minsize The minimum size of elevated read coverage that ProActive searches for. Default is 10000 base pairs.
#' @param maxsize The minimum size of elevated read coverage that ProActive searches for. Default is the close to the length of the contig being assessed.
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
shape_matcher <- function (microbialread_dataset, windowsize, minsize, maxsize, mode) {
  refnames <- unique(microbialread_dataset[,1])
  best_match_list <- list()
  filteredout_contigs <- rep(NA, length(refnames))
  reason <- rep(NA, length(refnames))
  A <- 1
  B <- 1
  C <- 1
  for (i in refnames) {
    if(B == floor(length(refnames)/4)){
      cat("A quarter of the way done with shape_matching \n")
    }
    if(B == floor(length(refnames)/2)){
      cat("Half of the way done with shape_matching \n")
    }
    if(B == floor((length(refnames)*3)/4)){
      cat("Almost done with shape_matching! \n")
    }
    B <- B+1
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == i),]
     if (microbial_subset[nrow(microbial_subset),3]< 30000) {
      filteredout_contigs[C] <- i
      reason[C] <- "Too Short"
      C <- C+1
      next
    } else if (microbial_subset[(order(microbial_subset[,2], decreasing=TRUE))[50],2] <= 10) {
      filteredout_contigs[C] <-  i
      reason[C] <- "Low read cov"
      C <- C+1
      next
    }
    microbial_subset <- windowsize_func(microbial_subset,windowsize, mode)
    no_transduction_best_match <- notransduction_func_WC(microbial_subset)
    prophage_off_left_best_match <- prophage_off_left_func_WC(microbial_subset, windowsize, minsize, maxsize)
    prophage_off_right_best_match <-  prophage_off_right_func_WC(microbial_subset, windowsize, minsize, maxsize)
    full_prophage_best_match <- full_prophage_func_WC(microbial_subset, windowsize, minsize, maxsize)
    Full_RCgap_best_match <- full_RCgap_func(microbial_subset, windowsize, minsize, maxsize)
    RCgap_offLeft_best_match <- RCgap_off_left_func(microbial_subset, windowsize, minsize, maxsize)
    RCgap_offRight_best_match <- RCgap_off_right_func(microbial_subset, windowsize, minsize, maxsize)
    best_match_summary <- list(no_transduction_best_match, prophage_off_left_best_match, prophage_off_right_best_match, full_prophage_best_match, Full_RCgap_best_match,RCgap_offLeft_best_match, RCgap_offRight_best_match)
    best_match_score_summary <- c(no_transduction_best_match[[1]],prophage_off_left_best_match[[1]], prophage_off_right_best_match[[1]], full_prophage_best_match[[1]], Full_RCgap_best_match[[1]], RCgap_offLeft_best_match[[1]], RCgap_offRight_best_match[[1]]) %>% as.numeric()
    best_match <- best_match_summary[[which(best_match_score_summary == min(best_match_score_summary))[1]]]
    best_match_list[[A]] <- c(best_match, i)
    A <- A+1
  }
  cov_check_predictions <- maxmin_coverage_check(best_match_list, microbialread_dataset, windowsize, mode)
  filteredout_contigsdf <- cbind.data.frame(filteredout_contigs, reason)
  filteredout_contigsdf <- na.omit(filteredout_contigsdf)
  shape_matching_summary <- list(cov_check_predictions[[1]], cov_check_predictions[[2]], cov_check_predictions[[3]], cov_check_predictions[[4]], filteredout_contigsdf)
  return(shape_matching_summary)
}

#add names to items in list


