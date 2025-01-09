#' 'chunk' long contigs
#'
#' Subset long contigs in metagenome pileup into chunks for pattern-matching
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' `chunkSize` determines the size (in bp) of each 'chunk'. Default is 100000.
#' @keywords internal
contigChunks <- function(pileup, chunkSize) {
  refNames <- unique(pileup[, 1])
  chunkedPileup <- data.frame(matrix(ncol = 3, nrow = nrow(pileup)))
  if (max(pileup$position) < chunkSize) {
    warning("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(pileup)
  }
  for (i in seq_along(refNames)){
    refName <- refNames[[i]]
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    if (pileupSubset[nrow(pileupSubset), 3] > chunkSize) {
      refNameChunk <- rep(NA, nrow(pileupSubset))
      chunkSizeDiv <- chunkSize / 100
      X <- 1
      Y <- chunkSizeDiv
      Z <- 1
      repeat {
        refNameChunk[c(X:Y)] <- rep(paste0(refName, "_chunk_", Z), chunkSizeDiv)
        if (Z == max(pileupSubset$position) %/% (chunkSize)) {
          break
        }
        X <- X + chunkSizeDiv
        Y <- Y + chunkSizeDiv
        Z <- Z + 1
      }
      remainIdx <- which(is.na(refNameChunk))[1]
      if (!is.na(remainIdx)) {
        contigChunk <- pileupSubset[remainIdx:nrow(pileupSubset), ]
        refNameChunk[c(remainIdx:nrow(pileupSubset))] <- rep(paste0(refName, "_chunk_", Z + 1), nrow(contigChunk))
      }
      chunkedPileupSubset <- cbind.data.frame(refNameChunk, pileupSubset$coverage, pileupSubset$position)
      colnames(chunkedPileupSubset) <- c("refName", "coverage", "position")
      NAIdx <- which(is.na(chunkedPileup[, 1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx + nrow(pileupSubset) - 1)), ] <- chunkedPileupSubset
    } else {
      NAIdx <- which(is.na(chunkedPileup[, 1]))[[1]]
      chunkedPileup[c(NAIdx:(NAIdx + nrow(pileupSubset) - 1)), ] <- pileupSubset
    }
  }
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
    warning("Your pileup is smaller than the chunkSize. The pileup will remain un-chunked.")
    return(pileup)
  }
  repeat {
    refName[c(X:Y)] <- rep(paste0(genomeRef, "_chunk_", Z), chunkSize)
    if (Z == max(pileup$position) %/% (chunkSize * 100)) break
    X <- X + chunkSize
    Y <- Y + chunkSize
    Z <- Z + 1
  }
  remainIdx <- which(is.na(refName))[1]
  genomeChunk <- pileup[remainIdx:nrow(pileup), ]
  refName[c(remainIdx:nrow(pileup))] <- rep(paste0(genomeRef, "_chunk_", Z + 1), nrow(genomeChunk))
  chunkedpileup <- cbind.data.frame(refName, pileup$coverage, pileup$position)
  colnames(chunkedpileup) <- c("refName", "coverage", "position")
  return(chunkedpileup)
}



#' Link pattern-matches on contig/genome chunks
#'
#' Detect partial gap/elevation pattern matches that fall on the edges of chunked
#' genomes/contigs that may be part of the same pattern prior to chunking
#'
#' @param bestMatchList A list containing pattern-match information associated with
#' all contigs/chunks classified by `ProActive()` pattern-matching
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param windowSize The number of basepairs to average read coverage values over.
#' @param mode Either "genome" or "metagenome"
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @importFrom stringr str_extract
#' @keywords internal
linkChunks <- function(bestMatchList, pileup, windowSize, mode, verbose){
  potLink <- rep(NA, length(bestMatchList))
  refNames <- vapply(seq_along(bestMatchList), function(i){bestMatchList[[i]][[8]]}, character(1))
  classifVector <- vapply(seq_along(bestMatchList), function(i){
    refName <- bestMatchList[[i]][[8]]
    classification <- bestMatchList[[i]][[7]]
    startPos <- bestMatchList[[i]][[4]]
    endPos <- bestMatchList[[i]][[5]]
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    if(grepl("chunk", refName) == FALSE) {
      "NoChunk"
    } else if(classification == "NoPattern") {
      "NoPat"
    } else if(startPos == 1) {
      "Left"
    } else if(endPos == nrow(pileupSubset)){
      "Right"
    } else {
      "Full"
    }
  }, character(1))
  classifDf <- cbind.data.frame(refNames, classifVector, potLink)
  rightIdxs <- which(classifDf[,2] == "Right")
  if(length(rightIdxs) != 0){
    for(p in seq_along(rightIdxs)){
      rightIdx <- rightIdxs[[p]]
      if(rightIdx == nrow(classifDf)) {classifDf[rightIdx,3] <- "no"}
      else if(classifDf[rightIdx + 1, 2] == "Left") {
        potLinkRef <- classifDf[rightIdx, 1]
        potLinkRef2 <- classifDf[rightIdx + 1, 1]
        sameChunk <- ifelse(gsub("_.*", "", potLinkRef) == gsub("_.*", "", potLinkRef2), "same", "notSame")
        seqChunk <- ifelse(as.numeric(str_extract(potLinkRef2, "[^_]+$")) == (as.numeric(str_extract(potLinkRef, "[^_]+$")) + 1) , "seq", "notSeq")
        classifDf[rightIdx,3] <- ifelse(sameChunk == "same" & seqChunk == "seq", "link", "no")
      } else {
        classifDf[rightIdx,3] <- "no"
      }
    }
    if(verbose == TRUE){
      if("link" %in% classifDf[,3]){
      linkIdxs <- which(classifDf[,3] == "link")
      lapply(seq_along(linkIdxs), function(x){
        linkIdx <- linkIdxs[[x]]
        message("Possible pattern-match continuity detected between ",
                                    classifDf[linkIdx,1], " and ", classifDf[linkIdx+1,1])
      })
    }
    }
  }
}
