#' Highly active/abundant prophage prediction function
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
#'shape_matching_results <- ProActive(whole_commreadcoverages, 1000)
#'}
ProActive <- function(microbial_readdataset, windowsize = 1000){
  start_time <- Sys.time()
  microbial_readdataset <- readcovdf_formatter(microbial_readdataset)
  print("Starting shape-matching")
  SM_predictions_summary <- shape_matcher_WC(microbial_readdataset, windowsize)
  SM_prophagepredictions_list <- SM_predictions_summary[[1]]
  SM_honmentions_list <- SM_predictions_summary[[2]]
  SM_none_predictions_list <- SM_predictions_summary[[3]]
  FullProphagePredictionList <- append(SM_prophagepredictions_list, SM_honmentions_list)
  filteredoutcontigs_df <- SM_predictions_summary[[4]]
  print("Identifying potential active prophages")
  Prediction_summary_df <- contig_prediction_summary_WC(FullProphagePredictionList)
  print("Determining sizes (bp) of potential active prophages")
  summary_table_prophagematchsize <- activeprophage_matchsize_checker(Prediction_summary_df, FullProphagePredictionList, windowsize)
  print("Finalizing output")
  final_summary_list<-list(summary_table_prophagematchsize, SM_prophagepredictions_list, SM_honmentions_list, filteredoutcontigs_df)
  end_time <- Sys.time()
  print(paste("Execuion time:", end_time-start_time))
  print(table(final_summary_list[[1]][,2]))
  print(paste(length(SM_prophagepredictions_list), "contigs with reliable predictions and", length(SM_honmentions_list), "contigs with unreliable predictions"))
  return(final_summary_list)
}

