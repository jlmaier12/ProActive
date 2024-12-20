#' Collects prophage classification shape-match information
#'
#' Collects shape information associated with all contigs predicted as containing prophages based on the results from the shape_matcher function.
#'
#' @param bestMatchList classifications made with shape_matcher function. classifications are stored as the first item in the bestMatchList.
#'
#' @keywords internal
allGapsAndElevs <- function(bestMatchList){
  A<-1
  gapElivList <- list()
  lapply(seq_along(bestMatchList), function(i){
    classification <-  bestMatchList[[i]][[7]]
    if(classification=="NoPattern") {
      return(NULL)
    }
    gapElivList[[A]] <<- bestMatchList[[i]]
    A <<- A+1
  })
  return(gapElivList)
}
