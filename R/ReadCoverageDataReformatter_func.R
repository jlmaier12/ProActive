#' Correctly formats input read coverage summary files.
#'
#' Places columns in correct order and renames columns. Clean the contig labels to remove excess informatio.
#'
#' @param read_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
readcovdf_formatter <- function(read_dataset, mode) {
  column_classes <- lapply(1:ncol(read_dataset), function(i) class(read_dataset[, i]))
  for (i in c(which(column_classes=="integer"))){
    if(mode=="genome"){
      if (read_dataset[1,i]==0){
        position_colindex <- i
      }
    } else {
    if (length(which(read_dataset[,i]==100))>1){
      position_colindex <- i
    }
    }
  }
  if (exists("position_colindex") == FALSE){
    cat("Error: Make sure you used binsize=100 when generating your pileup file.")
  }
  reformatted_readdataset <- cbind.data.frame(read_dataset[,which(column_classes == "character")],read_dataset[,which(column_classes == "numeric")], read_dataset[,position_colindex])
  colnames(reformatted_readdataset) <- c("ref_name", "coverage", "position")
  reformatted_readdataset$ref_name <- gsub("\\s.*", "", reformatted_readdataset$ref_name)
  reformatted_readdataset$ref_name <- str_extract(reformatted_readdataset$ref_name,regex("([^_]*_*[^_]*)"))
  return(reformatted_readdataset)
}

