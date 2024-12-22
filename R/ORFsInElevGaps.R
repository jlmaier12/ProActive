#' Detect ORFs in elevations and gaps
#'
#' Extracts subsets of the gffTSV associated with ORFs that fall within regions
#' of detected gapped or elevated read coverage.
#'
#' @param elevGapSummList A list containing pattern-match information associated with all
#' elevation and gap classifications. (i.e. no NoPattern classifications)
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param gffTSV Optional, a .gff file (TSV) containing ORFs associated with the .fasta
#' file used to generate the pileup.
#' @param mode Either "genome" or "metagenome"
#' @param chunkContigs TRUE or FALSE, If TRUE and `mode`="metagenome", contigs longer
#' than the `chunkSize` will be 'chunked' into smaller subsets and pattern-matching
#' will be performed on each subset. Default is FALSE.
#' @importFrom stringr str_detect str_extract_all regex
#' @importFrom dplyr bind_rows
#' @keywords internal
ORFsInElevGaps <- function(elevGapSummList, windowSize, gffTSV, mode, chunkContigs) {
  ORFlist <- list()
  colnames(gffTSV) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  if (TRUE %in% (str_detect(gffTSV[,9], regex('product', ignore_case = T)))){
    product <- str_extract_all(gffTSV [,9], regex("(?<=product=)[\\s\\S]*",ignore_case = T))
    gffTSV$geneproduct <- product
  }
  lapply(seq_along(elevGapSummList), function(i) {
    trueRefName <- elevGapSummList[[i]][[8]]
    if(mode == "metagenome"){
      refName <- elevGapSummList[[i]][[8]]
    if(chunkContigs==TRUE){
      refName <- gsub("_chunk_.*", "", refName)
    }
    gffTSV <- gffTSV[which(gffTSV[, 1] == refName), ]
    }
    startPos <- elevGapSummList[[i]][[4]] * windowSize
    endPos <- elevGapSummList[[i]][[5]] * windowSize
    matchRegion <- seq(startPos, endPos, 1)
    ORFs <- gffTSV[which(gffTSV[, 5] %in% matchRegion), ]
    ORFs$seqid <- rep(trueRefName, nrow(ORFs))
    if (nrow(ORFs) == 0) {
      return(NULL)
    }
    ORFs$Classification <- elevGapSummList[[i]][[7]]
    ORFlist[[i]] <<- ORFs
  })
  ORFSummTable <- bind_rows(ORFlist)
  return(ORFSummTable)
}
