#' shape-match size calculator
#'
#' Calculate the size, in base pairs, of the matching region for shapes detected in whole-community read coverages
#'
#' @param elevGapSummList A list containing shape information associated with all contigs containing a potential active prophages. Generated with the allprophages_func
#' @param windowSize The window size used to re-average read coverage datasets
#' @param ORFTable A gff file of ORFs
#' @keywords internal
ORFsInElevGaps <- function(elevGapSummList, windowSize, ORFTable){
  ORFlist <- list()
  colnames(ORFTable) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  lapply(seq_along(elevGapSummList),function(i){
    refName <- elevGapSummList[[i]][[8]]
    ORFSubset <- ORFTable[which(ORFTable[,1]==refName),]
    startPos <- elevGapSummList[[i]][[4]] * windowSize
    endPos <- elevGapSummList[[i]][[5]] * windowSize
    matchRegion <- seq(startPos, endPos, 1)
    ORFs <- ORFSubset[which(ORFSubset[,5] %in% matchRegion),]
    if(nrow(ORFs)==0){
      return(NULL)
    }
    ORFs$Classification <- elevGapSummList[[i]][[7]]
    ORFlist[[i]] <<- ORFs
  })
  ORFSummTable <- bind_rows(ORFlist)
  return(ORFSummTable)
}
