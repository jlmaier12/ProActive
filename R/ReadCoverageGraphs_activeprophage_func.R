#' Plot read coverage graphs of contigs with predicted active prophages in the whole-community read coverages
#'
#' Plot the read coverages of a contig and its associated shape-match for each contig with a predicted active prophage.
#'
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param prophagepredictions A list containing shape information associated with all contigs containing a potential active prophage. The third item in the list output from the main prophage_activity_finder_func
#' @param final_summary_table Prediction summary table with active prophage information. The first item in the list output from the main prophage_activity_finder_func function.
#' @param windowsize The window size used to re-average read coverage datasets
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ProActivePredictionPlots(whole_commreadcovs, shapematch_result[[3]],shapematch_result[[1]],1000)
#' }
ProActivePredictionPlots <- function(microbialread_dataset, prophagepredictions,final_summary_table, windowsize = 1000) {
  position <- coverage <- NULL
  windowsize <- windowsize
  microbialread_dataset <- readcovdf_formatter(microbialread_dataset)
  for (i in seq(1,length(prophagepredictions),1)) {
    ref_name <- prophagepredictions[[i]][[7]]
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == ref_name),]
    microbial_subset <- windowsize_func(microbial_subset,windowsize)
    microbial_subset[is.nan.data.frame(microbial_subset)] <- 0
    match_info <- final_summary_table[which(final_summary_table[,1]==ref_name),]
    shape <- shape_builder_func_WC(microbial_subset, prophagepredictions, i)
    shape_match <- cbind(microbial_subset, shape)
    match_length <- match_info[,4]
    print(ggplot(data=shape_match, aes(x=position, y=coverage))+
            geom_area(fill="deepskyblue3") +
            geom_line(y=shape, size=1)+
            labs(title=ref_name,subtitle=paste("Matching-region size (bp):", match_length), x=" ") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15)))
  }
}
