#' Search for gene annotations on classified contigs/chunks
#'
#' Search contigs classified with ProActive for gene-annotations that match a provided
#' key-word(s). Outputs read coverage plots for contigs/chunks with matching annotations.
#'
#' @param ProActiveResults The output from `ProActive()`.
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param inGapOrElev TRUE or FALSE. If TRUE, only search for gene-annotations in
#' the gap/elevation region of the pattern-match. Default is FALSE (i.e search the
#' entire contig/chunk for the gene annotation key-words)
#' @param bpRange If `inGapOrElev` = TRUE, the user may specify the region (in base pairs) that should
#' be searched to the left and right of the gap/elevation region. Default is 0.
#' @param gffTSV A .gff file (TSV) containing gene predictions associated with the .fasta
#' file used to generate the pileup.
#' @param geneOrProduct "gene" or "product". Search for keyWords associated with genes or gene products.
#' @param keyWords The keyWord(s) to search for. Case independent. Searches will return the string
#' that contains the matching keyWord. KeyWord(s) must be in quotes, comma-separated, and surrounded by
#' c() i.e( c("antibiotic", "resistance", "drug") )
#' @param elevFilter Optional, only plot results with pattern-matches that achieved an
#' elevation ratio (max/min) greater than the specified values. Default is no filter.
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#' output to. A folder will be made within the provided directory to store
#' results.
#' @returns list of ggplot objects
#' @importFrom dplyr %>%
#' @export
#' @examples
#' geneAnnotMatches <- geneAnnotationSearch(sampleMetagenomeResults, sampleMetagenomePileup,
#'                                           sampleMetagenomegffTSV, geneOrProduct="product",
#'                                           keyWords=c("toxin", "drug", "resistance", "phage"))
geneAnnotationSearch <- function(ProActiveResults, pileup, gffTSV,
                                 geneOrProduct, keyWords, inGapOrElev = FALSE,
                                 bpRange = 0, elevFilter, saveFilesTo, verbose = TRUE) {
  if(bpRange != 0 & inGapOrElev == FALSE) {stop("Cannot set bpRange if inGapOrElev = FALSE")}
  if(verbose){message("Cleaning gff file...")}
  colnames(gffTSV) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  if (tolower(geneOrProduct) == "product") {
    if (TRUE %in% (str_detect(gffTSV$attributes, regex('product', ignore_case = TRUE)))){
    product <- str_extract(gffTSV$attributes, regex("(?<=product=)[\\s\\S]*", ignore_case = TRUE))
    gffTSV$product <- product
    } else { stop("No products detected in attributes column of gff file.")}
  } else {
    if (TRUE %in% (str_detect(gffTSV$attributes, regex('gene', ignore_case = TRUE)))){
    gene <- str_extract(gffTSV$attributes, regex("(?<=gene=)[\\s\\S]*", ignore_case = TRUE)) %>% str_extract(".+?(?=;)")
    gffTSV$gene <- gene
    } else { stop("No genes detected in attributes column of gff file.")}
  }
  elevFilter <- ifelse(missing(elevFilter), NA, elevFilter)
  summaryTable <- ProActiveResults[[1]]
  patternMatches <- ProActiveResults[[3]]
  mode <- ProActiveResults[[5]][[2]]
  chunkSize <- ProActiveResults[[5]][[3]]
  contigChunk <- ProActiveResults[[5]][[4]]
  windowSize <- ProActiveResults[[5]][[1]]
  if(verbose){message("Cleaning pileup file...")}
  pileup <- pileupFormatter(pileup, mode)
  if (mode == "genome") {
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if (mode == "metagenome" & contigChunk) {
    pileup <- contigChunks(pileup, chunkSize)
  }
  if(verbose){message("Searching for matching annotations...")}
  plots <- lapply(seq_along(patternMatches), function(i){
                  trueRefName <- patternMatches[[i]][[8]]
                  classification <- patternMatches[[i]][[7]]
                  elevRatio <- ifelse(classification == "NoPattern", 0, patternMatches[[i]][[6]])
                  if (is.na(elevFilter) == FALSE & elevRatio < elevFilter) {return(NULL)}
                  if (inGapOrElev & classification == "NoPattern") {return(NULL)}
                  pileupSubset <- pileup[which(pileup[,1] == trueRefName),]
                  pileupRegion <- seq(pileupSubset[1, 3], pileupSubset[nrow(pileupSubset), 3], 1)
                  if(mode == "metagenome"){
                    refName <- patternMatches[[i]][[8]]
                    if(contigChunk){
                    refName <- gsub("_chunk_.*", "", refName)
                    }
                    gffTSV <- gffTSV[which(gffTSV[, 1] == refName), ]
                  }
                  geneAnnotSubset <- gffTSV[which(gffTSV$end %in% pileupRegion), ]
                  colIdx <- ifelse(tolower(geneOrProduct) == "gene",
                                     which(colnames(gffTSV) == "gene"),
                                     which(colnames(gffTSV) == "product"))
                  if (TRUE %in% (str_detect(geneAnnotSubset[, colIdx], regex(paste(keyWords, collapse="|"), ignore_case = TRUE)))) {
                    startbpRange <-  endbpRange <- NULL
                    if (inGapOrElev){
                       if(mode == "genome"){
                         chunkNumber <- as.numeric(str_extract(trueRefName, "(?<=\\_)\\d+$")) - 1
                         startPos <- (patternMatches[[i]][[4]] * windowSize) + (chunkNumber * chunkSize)
                         endPos <- (patternMatches[[i]][[5]] * windowSize) + (chunkNumber * chunkSize)
                       } else {
                        startPos <- patternMatches[[i]][[4]] * windowSize
                        endPos <- patternMatches[[i]][[5]] * windowSize
                      }
                      endbpRange <- ifelse ((endPos + bpRange > (pileupSubset[nrow(pileupSubset), 3])),
                                            (pileupSubset[nrow(pileupSubset), 3]), (endPos + bpRange))
                      startbpRange <- ifelse((startPos - bpRange < pileupSubset[1, 3] ), pileupSubset[1, 3], (startPos - bpRange))
                      matchRegion <- seq(startbpRange, endbpRange, 1)
                      geneAnnotSubset <- gffTSV[which(gffTSV$end %in% matchRegion), ]
                      if (!(TRUE %in% (str_detect(geneAnnotSubset[, colIdx], regex(paste(keyWords, collapse="|"), ignore_case = TRUE))))) {
                        return(NULL)
                      }
                    }
                    plot <- geneAnnotationPlot(geneAnnotSubset, keyWords,
                                                 pileupSubset, colIdx, startbpRange, endbpRange,
                                                 elevRatio, patternMatches[[i]], windowSize, chunkSize, mode)
                  } else {return(NULL)}
    plot
  })
  plots <- (plots[!vapply(plots, is.null, logical(1))])
  refNames <- vapply(seq_along(plots), function(i){
    plot <- plots[[i]]
    plotdata <- ggplot_build(plot)
    contig <- plotdata$plot$labels$title %>% str_extract(".+?(?=[:space:])")
    contig
  }, character(1))
  names(plots) <- refNames
  if(verbose){message(length(plots), " contigs/chunks have gene annotations that match one or more of the provided keyWords")}
  if (missing(saveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\ProActiveGeneAnnotMatchPlots")),
           dir.create(paste0(saveFilesTo, "\\ProActiveGeneAnnotMatchPlots")),
           stop(
             "'ProActiveGeneAnnotMatchPlots' already exists in the provided directory"
           )
    )
    lapply(
      names(plots),
      function(X) {
        ggsave(
          filename = paste0(
            saveFilesTo,
            "\\ProActiveGeneAnnotMatchPlots\\", X, ".png"
          ),
          plot = plots[[X]],
          width = 8,
          height = 4
        )
      }
    )
    return(plots)
  } else {
    return(plots)
  }
}
