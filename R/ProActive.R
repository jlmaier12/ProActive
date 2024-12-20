#' Detect elevations and gaps in mapped read coverage patterns.
#'
#' Performs read coverage pattern-matching and summarizes the results into a list.
#' The first item in the list is a table consisting of the summary information of all the
#' contigs that passed through shape-matching (i.e were not filtered out). The second item
#' in the list is a table consisting of the summary information of all contigs that were
#' predicted as containing a potential active prophage. The third item in the list contains
#' the best shape-match information associated with each contig in the previous table. The
#' fourth and final object in the list is a table containing the contigs that were filtered
#' out prior to shape_matching and the reason why.
#'
#'@param pileup A pileup file for either a genome or metagenome containing mapped read coverages averaged over 100 bp window sizes
#'@param mode Either "genome" or "metagenome"
#'@param gffTSV Optional,
#'@param windowSize The number of base pairs to average coverage values over. Must be either 200, 500, or 1000. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000 is the default.
#'@param chunkContigs
#'@param minSize The minimum size (in bp) of elevated/gapped read coverage to detect. Default is 10000.
#'@param maxSize The minimum size (in bp) of elevated/gapped read coverage to detect. Default is NA (i.e no maximum).
#'@param minContigLength
#'@param chunkSize If mode="genome", ProActive will split your genome into chunks to search for elevated read coverage. Since ProActive can only identify one elevated read coverage pattern per 'chunk', choosing smaller chunks will increase sensitivity AND processing time. Default is 100000 (bp per chunk).
#'@param IncludeNoPatterns TRUE or FALSE. Default is FALSE. TRUE if you would like to include the contigs classified as 'none' in your final pattern lists. This can be useful if you are comparing ProActive results between conditions.
#'@param saveFilesTo Optional, Provide a path to the directory you wish to save output to. A folder will be made within the provided directory to store results.
#'@export
#'
#'@examples
#' \dontrun{
#' ## Metagenome mode
#' data(sampleMetagenomePileup)
#' data(sampleMetagenomegffTSV)
#'ProActive_results <- ProActive(pileup=sampleMetagenomePileup,
#'                               mode="metagenome",
#'                               gffTSV=sampleMetagenomegffTSV)
#'
#' ## Genome mode
#'
#'}
ProActive <- function(pileup, mode, gffTSV, windowSize = 1000, chunkContigs= FALSE, minSize=10000, maxSize=Inf, minContigLength=30000, chunkSize=100000, IncludeNoPatterns=FALSE, saveFilesTo){
  ##error catching
  if ((chunkSize %% 100) >0) {
    stop("chunkSize must be divisible by 100")
  }
  if (!(windowSize %in% list(100, 200, 500, 1000))) {
    stop("windowSize must be either 100, 200, 500, or 1000 bp!")
  }
  if (minContigLength <= 25000) {
    stop("minContigLength must be at least 25,000 bp for pattern-matching!")
  }
  if (minSize <= 1000) {
    stop("minBlockSize must be greater than 1000 bp!")
  }
  if (abs(pileup[1, 3] - pileup[2, 3]) != 100)  {
    stop("Pileup file MUST have a windowSize/binsize of 100!")
  }
  startTime <- Sys.time()
  message("Preparing input file for pattern-matching...")
    pileup <- pileupFormatter(pileup, mode)
  if(mode=="genome"){
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if(chunkContigs == TRUE){
    pileup <- contigChunks(pileup, chunkSize)
  }
  message("Starting pattern-matching...")
  patternMatchSummary <- patternMatcher(pileup, windowSize, minSize, maxSize, mode, minContigLength)
  if (IncludeNoPatterns==TRUE) {
    classifList <- patternMatchSummary[[1]]
  } else {
    classifList <- removeNoPatterns(patternMatchSummary[[1]])
  }
  filteredOutContigsDf <- patternMatchSummary[[2]]
  message("Summarizing pattern-matching results")
  summaryTable <- classifSumm(pileup, patternMatchSummary[[1]], windowSize, mode)
  if (missing(gffTSV) == FALSE) {
    message("Finding ORFs in elevated or gapped regions of read coverage...")
    elevGapSummList <- removeNoPatterns(patternMatchSummary[[1]])
    ORFSummTable <- ORFsInElevGaps(elevGapSummList, windowSize, gffTSV)
  }
  message("Finalizing output")
  endTime <- Sys.time()
  duration <- difftime(endTime, startTime)
  message("Execution time: ", round(duration[[1]], 2), units(duration))
  message(
    length(which(
      filteredOutContigsDf[, 2] == "Low read cov"
    )),
    " contigs were filtered out based on low read coverage"
  )
  message(
    length(which(
      filteredOutContigsDf[, 2] == "Too Short"
    )),
    " contigs were filtered out based on length (< minContigLength)"
  )
    arguments <- list(windowSize, mode, chunkSize, chunkContigs)
    finalSummaryList<-list(summaryTable, classifList, filteredOutContigsDf, arguments)
    names(finalSummaryList) <- c("SummaryTable", "PatternMatches","FilteredOut", "Arguments")
    if(missing(gffTSV) == FALSE){
      finalSummaryList <- c(finalSummaryList, list(ORFSummTable))
      names(finalSummaryList)[5] <- "ORFTable"
    }
    table <- (table(summaryTable[, 2]))
    message(paste0(capture.output(table), collapse = "\n"))
  if (missing(saveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(SaveFilesTo, "\\ProActiveOutput")),
           dir.create(paste0(SaveFilesTo, "\\ProActiveOutput")),
           stop(
             "'ProActiveOutput' folder exists already
                in the provided directory"
           )
    )
    if(missing(gffTSV) == FALSE){
      write.table(
        ORFSummTable,
        file = paste0(
          SaveFilesTo,
          "\\ProActiveOutput\\ProActiveORFResults.csv"
        ),
        sep = ",",
        row.names = FALSE
      )
    }
    write.table(
      summaryTable,
      file = paste0(
        SaveFilesTo,
        "\\ProActiveOutput\\ProActiveSummaryTable.csv"
      ),
      sep = ",",
      row.names = FALSE
    )
    write.table(
      filteredOutContigsDf,
      file = paste0(
        SaveFilesTo,
        "\\ProActiveOutput\\ProActiveFilteredOutContigs.csv"
      ),
      sep = ",",
      row.names = FALSE
    )
    return(finalSummaryList)
  } else {
    return(finalSummaryList)
  }
}



