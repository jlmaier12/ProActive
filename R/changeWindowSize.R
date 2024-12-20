#' Change the read coverage window size
#'
#' Re-averages window sizes of read coverage averages. Start with 100bp windows always. Cannot make window size less than 100bp.
#'
#' @param pileup A read coverage dataset that has been cleaned and reformatted by the readcovdf_formatter function
#' @param windowSize The number of base pairs to average coverage values over
#' @param mode Either "genome" or "metagenome"
#'
#' @keywords internal
changewindowSize <- function(pileup, windowSize, mode){ #anything divisible by 100, pileup should be a single contig
  coverage <- vector()
  X <- 0
  Y <- windowSize/100
  repeat{
    coverage <- c(coverage, mean(pileup[c(X:Y),2]))
    X <- X+(windowSize/100)
    Y <- Y+(windowSize/100)
    if (Y > nrow(pileup)) break
  }
  if(mode=="genome"){
    position <- seq(pileup[1,3], pileup[nrow(pileup),3], length.out=length(coverage))
  } else {
  position <- seq(windowSize, length(coverage)*windowSize, windowSize)
  }
  refName <- rep(pileup[1,1], length(position))
  newdataset <- cbind.data.frame(refName, coverage, position) %>% as.data.frame()
  newdataset[do.call(cbind, lapply(newdataset, is.nan))] <- 0
  return(newdataset)
}
