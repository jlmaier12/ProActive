#' Main transduction prediction function
#'
#' The main function that ties together all the other functions. This function performs all the shape-matching and summarizes the results into a list. The first item in the list is a table consisting of the summary information of all the contigs that passed through shape-matching (i.e were not filtered out). The second item in the list is a table consisting of the summary information of all contigs that were predicted as containing a potential active prophage. The third item in the list contains the best shape-match information associated with each contig in the previous table. The fourth and final object in the list is a table containing the contigs that were filtered out prior to shape_matching and the reason why.
#'
#'@param microbial_readdataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#'@param windowsize The number of base pairs to average coverage values over. For compatibility with shape-matching, window sizes can only be 200, 500, 1000 and 2000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowsize is the default.
#'
#'@export
#'
#'@examples
#' \dontrun{
#'shape_matching_results <- prophage_activity_finder_func(whole_commreadcoverages, 1000)
#'}
prophage_activity_finder_func <- function(microbial_readdataset, windowsize = 1000){
  start_time <- Sys.time()
  windowsize <- windowsize
  microbial_readdataset <- readcovdf_formatter(microbial_readdataset)
  print("Starting shape-matching")
  SM_predictions_summary <- shape_matcher_WC(microbial_readdataset, windowsize)
  SM_predictions <- SM_predictions_summary[[1]]
  filteredoutcontigs <- SM_predictions_summary[[2]]
  print("Identifying potential active prophages")
  prophage_predictions_list <- allprophages_func_WC(SM_predictions)
  Prediction_summary_df <- contig_prediction_summary_WC(SM_predictions)
  print("Determining sizes (bp) of potential active prophages")
  summary_table_matchsize <- activeprophage_matchsize_checker(Prediction_summary_df, prophage_predictions_list, windowsize)
  print("Finalizing output")
  cleaned_summary_table <- summary_table_matchsize[which(summary_table_matchsize[,2]=="Prophage"),]
  final_summary_list<-list(summary_table_matchsize,cleaned_summary_table, prophage_predictions_list, filteredoutcontigs)
  end_time <- Sys.time()
  print(paste("Execuion time:", end_time-start_time))
  print(table(final_summary_list[[1]][,2]))
  return(final_summary_list)
}

