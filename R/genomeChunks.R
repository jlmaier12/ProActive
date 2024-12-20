#' Break genome pileup into chunks for pattern-matching
#'
#' @param genomePileup Pileup file for read coverages mapped to a genome
#' @param chunkSize Value to chunk genome into. Should be divisible by 100. Default is 100,000bp per chunk.
#'
#' @keywords internal
genomeChunks <- function(genomePileup, chunkSize) {
  refName <- rep(NA, nrow(genomePileup))
  genomeRef <- genomePileup[1,1]
  chunkSize <- chunkSize/100
  X <- 1
  Y <- chunkSize
  Z <- 1
  if(max(genomePileup$position) < chunkSize){
    message("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(genomePileup)
  }
  repeat {
  refName[c(X:Y)] <- rep(paste0(genomeRef,"_chunk_", Z), chunkSize)
  if(Z==max(genomePileup$position)%/%(chunkSize*100)) break
  X <- X+chunkSize
  Y <- Y+chunkSize
  Z <- Z+1
  }
  remain_index <- which(is.na(refName))[1]
  genomeChunk <- genomePileup[remain_index:nrow(genomePileup),]
  refName[c(remain_index:nrow(genomePileup))] <- rep(paste0(genomeRef, "_chunk_", Z+1), nrow(genomeChunk))
  chunkedGenomePileup <- cbind.data.frame(refName, genomePileup$coverage, genomePileup$position)
  colnames(chunkedGenomePileup) <- c("refName", "coverage", "position")
  return(chunkedGenomePileup)
}


