#' Correctly formats pileup files.
#'
#' Places columns in correct order and renames columns. Cleans the contig labels
#' to remove excess information after whitespace.
#'
#' @param pileup A table containing contig names, read coverages averaged over
#'   100 bp windows,and contig positions
#' @param mode Either "genome" or "metagenome"
#' @return dataframe
#' @keywords internal
pileupFormatter <- function(pileup, mode) {
  colClasses <-
    lapply(seq_along(pileup), function(i) {
      class(pileup[, i])
    })
  for (i in c(which(colClasses == "integer"))) {
    if(mode=="genome"){
      if (pileup[1,i]==0){
        posColIdx <- i
      }
  } else {
    if (length(which(pileup[, i] == 100)) > 1) {
      posColIdx <<- i
    }
  }
  }
  cleanPileup <-
    cbind.data.frame(
      pileup[, which(colClasses == "character")],
      pileup[, which(colClasses == "numeric")],
      pileup[, posColIdx]
    )
  colnames(cleanPileup) <- c("refName", "coverage", "position")
  cleanPileup$refName <- gsub("\\s.*", "", cleanPileup$refName)
  return(cleanPileup)
}

