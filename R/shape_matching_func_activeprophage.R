#' Shape-matching function for active prophages in whole-community read coverage
#'
#' Creates the microbial_subset, representative of one contig, that is used as input for each individual shape-building function. After the information associated with the best match for each shape is obtained, the shape with the lowest mean absolute difference (score) is chosen as the prediction for the contig being assessed.
#'
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets
shape_matcher_WC <- function (microbialread_dataset, windowsize) {
  windowsize <- windowsize
  microbialread_dataset <- readcovdf_formatter(microbialread_dataset)
  refnames <- unique(microbialread_dataset[,1])
  best_match_list <- list()
  filteredout_contigs <- rep(NA, length(refnames))
  A <- 1
  B <- 1
  C <- 1
  for (i in refnames) {
    if(B == floor(length(refnames)/4)){
      print("A quarter of the way done with shape_matching")
    }
    if(B == floor(length(refnames)/2)){
      print("Half of the way done with shape_matching")
    }
    if(B == floor((length(refnames)*3)/4)){
      print("Almost done with shape_matching!")
    }
    B <- B+1
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == i),]
    if (microbial_subset[nrow(microbial_subset),3]< 30000) {
      filteredout_contigs[C] <- i
      C <- C+1
      next
    } else if (microbial_subset[(order(microbial_subset[,2], decreasing=TRUE))[50],2] <= 10) {
      filteredout_contigs[C] <-  i
      C <- C+1
      next
    }
    microbial_subset <- windowsize_func(microbial_subset,windowsize)
    microbial_subset[is.nan.data.frame(microbial_subset)] <- 0
    no_transduction_best_match <- notransduction_func_WC(microbial_subset)
    prophage_off_left_best_match <- prophage_off_left_func_WC(microbial_subset, windowsize)
    prophage_off_right_best_match <-  prophage_off_right_func_WC(microbial_subset, windowsize)
    full_prophage_best_match <- full_prophage_func_WC(microbial_subset, windowsize)
    best_match_summary <- list(no_transduction_best_match, prophage_off_left_best_match, prophage_off_right_best_match, full_prophage_best_match)
    best_match_score_summary <- c(no_transduction_best_match[[1]],prophage_off_left_best_match[[1]], prophage_off_right_best_match[[1]], full_prophage_best_match[[1]]) %>% as.numeric()
    best_match <- best_match_summary[[which(best_match_score_summary == min(best_match_score_summary))[1]]]
    best_match_list[[A]] <- append(best_match, i)
    A <- A+1
  }
  filteredout_contigs <- filteredout_contigs[!is.na(filteredout_contigs)]
  shape_matching_summary <- list(best_match_list, filteredout_contigs)
  return(shape_matching_summary)
}


