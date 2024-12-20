#' Break long contigs in metagenome pileup into chunks for pattern-matching
#'
#' @param pileup Pileup file for read coverages mapped to a genome
#' @param chunkSize Value to chunk genome into. Should be divisible by 100. Default is 100,000bp per chunk.
#'
#' @keywords internal
contigChunks <- function(pileup, chunkSize) {
  refNames <- unique(pileup[,1])
  chunkedPileup = data.frame(matrix(ncol = 3, nrow = nrow(pileup)))
  if(max(pileup$position) < chunkSize){
    message("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(pileup)
  }
  lapply(seq_along(refNames), function(i){
    refName <<- refNames[[i]]
    pileupSubset <<- pileup[which(pileup[,1]==refName),]
    if (pileupSubset[nrow(pileupSubset),3] > chunkSize) {
      refNameChunk <<- rep(NA, nrow(pileupSubset))
      chunkSizeDiv <<- chunkSize/100
      X <<- 1
      Y <<- chunkSizeDiv
      Z <<- 1
      repeat {
        refNameChunk[c(X:Y)] <<- rep(paste0(refName,"_chunk_", Z), chunkSizeDiv)
        if(Z==max(pileupSubset$position)%/%(chunkSize)) {
          break
        }
        X <<- X+chunkSizeDiv
        Y <<- Y+chunkSizeDiv
        Z <<- Z+1
      }
      remainIdx <- which(is.na(refNameChunk))[1]
      if(!is.na(remainIdx)){
      contigChunk <- pileupSubset[remainIdx:nrow(pileupSubset),]
      refNameChunk[c(remainIdx:nrow(pileupSubset))] <- rep(paste0(refName,"_chunk_", Z+1), nrow(contigChunk))
      }
      chunkedPileupSubset <- cbind.data.frame(refNameChunk, pileupSubset$coverage, pileupSubset$position)
      colnames(chunkedPileupSubset) <- c("refName", "coverage", "position")
      NAIdx <- which(is.na(chunkedPileup[,1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx+nrow(pileupSubset)-1)),] <<- chunkedPileupSubset
    } else {
      NAIdx <- which(is.na(chunkedPileup[,1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx+nrow(pileupSubset)-1)),] <<- pileupSubset
    }
  })
  return(chunkedPileup)
}
