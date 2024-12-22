#' 'chunk' long contigs
#'
#' Subset long contigs in metagenome pileup into chunks for pattern-matching
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' `chunkSize` determines the size (in bp) of each 'chunk'. Default is 50000.
#' @keywords internal
contigChunks <- function(pileup, chunkSize) {
  refNames <- unique(pileup[, 1])
  chunkedPileup <- data.frame(matrix(ncol = 3, nrow = nrow(pileup)))
  if (max(pileup$position) < chunkSize) {
    message("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(pileup)
  }
  lapply(seq_along(refNames), function(i) {
    refName <<- refNames[[i]]
    pileupSubset <<- pileup[which(pileup[, 1] == refName), ]
    if (pileupSubset[nrow(pileupSubset), 3] > chunkSize) {
      refNameChunk <<- rep(NA, nrow(pileupSubset))
      chunkSizeDiv <<- chunkSize / 100
      X <<- 1
      Y <<- chunkSizeDiv
      Z <<- 1
      repeat {
        refNameChunk[c(X:Y)] <<- rep(paste0(refName, "_chunk_", Z), chunkSizeDiv)
        if (Z == max(pileupSubset$position) %/% (chunkSize)) {
          break
        }
        X <<- X + chunkSizeDiv
        Y <<- Y + chunkSizeDiv
        Z <<- Z + 1
      }
      remainIdx <- which(is.na(refNameChunk))[1]
      if (!is.na(remainIdx)) {
        contigChunk <- pileupSubset[remainIdx:nrow(pileupSubset), ]
        refNameChunk[c(remainIdx:nrow(pileupSubset))] <- rep(paste0(refName, "_chunk_", Z + 1), nrow(contigChunk))
      }
      chunkedPileupSubset <- cbind.data.frame(refNameChunk, pileupSubset$coverage, pileupSubset$position)
      colnames(chunkedPileupSubset) <- c("refName", "coverage", "position")
      NAIdx <- which(is.na(chunkedPileup[, 1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx + nrow(pileupSubset) - 1)), ] <<- chunkedPileupSubset
    } else {
      NAIdx <- which(is.na(chunkedPileup[, 1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx + nrow(pileupSubset) - 1)), ] <<- pileupSubset
    }
  })
  X <- NULL
  Y <- NULL
  Z <- NULL
  refName <- NULL
  pileupSubset <- NULL
  chunkSizeDiv <- NULL
  return(chunkedPileup)
}

#' 'chunk' genomes
#'
#' Subset genome pileup into chunks for pattern-matching
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' `chunkSize` determines the size (in bp) of each 'chunk'. Default is 50000.
#' @keywords internal
genomeChunks <- function(pileup, chunkSize) {
  refName <- rep(NA, nrow(pileup))
  genomeRef <- pileup[1, 1]
  chunkSize <- chunkSize / 100
  X <- 1
  Y <- chunkSize
  Z <- 1
  if (max(pileup$position) < chunkSize) {
    message("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(pileup)
  }
  repeat {
    refName[c(X:Y)] <- rep(paste0(genomeRef, "_chunk_", Z), chunkSize)
    if (Z == max(pileup$position) %/% (chunkSize * 100)) break
    X <- X + chunkSize
    Y <- Y + chunkSize
    Z <- Z + 1
  }
  remain_index <- which(is.na(refName))[1]
  genomeChunk <- pileup[remain_index:nrow(pileup), ]
  refName[c(remain_index:nrow(pileup))] <- rep(paste0(genomeRef, "_chunk_", Z + 1), nrow(genomeChunk))
  chunkedpileup <- cbind.data.frame(refName, pileup$coverage, pileup$position)
  colnames(chunkedpileup) <- c("refName", "coverage", "position")
  return(chunkedpileup)
}
