#' Detect elevations and gaps in mapped read coverage patterns.
#'
#' Performs read coverage pattern-matching and summarizes the results into a list.
#'  The first list item summarizes the pattern-matching results. The second list
#' item is the 'cleaned' version of the summary table with all the 'noPattern'
#' classifications removed. (i.e were not filtered out). The third list item contains
#' the pattern-match information needed for pattern-match visualization with
#' `plotProActiveResults()`. The fourth list item is a table containing all the contigs that
#' were filtered out prior to pattern-matching. The fifth list item contains arguments used during
#' pattern-matching (windowSize, mode, chunkSize, chunkContigs). If the user provides a
#' gffTSV files, then the last list is a table consisting of ORFs found within
#' the detected gaps and elevations in read coverage.
#'
#' @param pileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param mode Either "genome" or "metagenome"
#' @param gffTSV Optional, a .gff file (TSV) containing gene predictions associated with the .fasta
#' file used to generate the pileup.
#' @param windowSize The number of basepairs to average read coverage values over.
#' Options are 100, 200, 500, 1000 ONLY. Default is 1000.
#' @param chunkContigs TRUE or FALSE, If TRUE and `mode`="metagenome", contigs longer
#' than the `chunkSize` will be 'chunked' into smaller subsets and pattern-matching
#' will be performed on each subset. Default is FALSE.
#' @param minSize The minimum size (in bp) of elevation or gap patterns. Default is 10000.
#' @param maxSize The maximum size (in bp) of elevation or gap patterns. Default is NA
#' (i.e. no maximum).
#' @param minContigLength The minimum contig/chunk size (in bp) to perform pattern-matching
#' on. Default is 25000.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' `chunkSize` determines the size (in bp) of each 'chunk'. Default is 50000.
#' @param IncludeNoPatterns TRUE or FALSE, If TRUE the noPattern pattern-matches will
#' be included in the ProActive PatternMatches output list. If you would like to visualize
#' the noPattern pattern-matches in `plotProActiveResults()`, this should be set to TRUE.
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#' output to. A folder will be made within the provided directory to store
#' results.
#' @importFrom utils capture.output write.table
#' @return A list containing 6 objects described in the function description.
#' @export
#' @examples
#' ## Metagenome mode with gffTSV
#' metagenome_results <- ProActive(
#'   pileup = sampleMetagenomePileup,
#'   mode = "metagenome",
#'   gffTSV = sampleMetagenomegffTSV
#' )
#'
#' ## Genome mode without gffTSV
#' genome_results <- ProActive(
#'   pileup = exampleGenomePileupSubset,
#'   mode = "genome"
#' )
#'
#' ##gffTSV is optional!
ProActive <- function(pileup, mode, gffTSV, windowSize = 1000, chunkContigs = FALSE,
                      minSize = 10000, maxSize = Inf, minContigLength = 30000,
                      chunkSize = 100000, IncludeNoPatterns = FALSE, saveFilesTo) {
  ## error catching
  if ((chunkSize %% 100) > 0) {
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
  if (abs(pileup[1, 3] - pileup[2, 3]) != 100) {
    stop("Pileup file MUST have a windowSize/binsize of 100!")
  }
  startTime <- Sys.time()
  message("Preparing input file for pattern-matching...")
  pileup <- pileupFormatter(pileup, mode)
  if (mode == "genome") {
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if (chunkContigs == TRUE) {
    pileup <- contigChunks(pileup, chunkSize)
  }
  message("Starting pattern-matching...")
  patternMatchSummary <- patternMatcher(pileup, windowSize, minSize, maxSize, mode, minContigLength)
  if (IncludeNoPatterns == TRUE) {
    classifList <- patternMatchSummary[[1]]
  } else {
    classifList <- removeNoPatterns(patternMatchSummary[[1]])
  }
  filteredOutContigsDf <- patternMatchSummary[[2]]
  message("Summarizing pattern-matching results")
  summaryTable <- classifSumm(pileup, patternMatchSummary[[1]], windowSize, mode)
  if (missing(gffTSV) == FALSE) {
    message("Finding gene predictions in elevated or gapped regions of read coverage...")
    elevGapSummList <- removeNoPatterns(patternMatchSummary[[1]])
    GPSummTable <- GPsInElevGaps(elevGapSummList, windowSize, gffTSV, mode, chunkContigs)
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
  cleanSummaryTable <- summaryTable[-which(summaryTable[,2]=="NoPattern"),]
  finalSummaryList <- list(summaryTable, cleanSummaryTable, classifList, filteredOutContigsDf, arguments)
  names(finalSummaryList) <- c("SummaryTable", "CleanSummaryTable", "PatternMatches", "FilteredOut", "Arguments")
  if (missing(gffTSV) == FALSE) {
    finalSummaryList <- c(finalSummaryList, list(GPSummTable))
    names(finalSummaryList)[6] <- "GenePredictTable"
  }
  table <- (table(summaryTable[, 2]))
  message(paste0(capture.output(table), collapse = "\n"))
  if (missing(saveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\ProActiveOutput")),
      dir.create(paste0(saveFilesTo, "\\ProActiveOutput")),
      stop(
        "'ProActiveOutput' exists already in the provided directory"
      )
    )
    if (missing(gffTSV) == FALSE) {
      GPSummTable <- apply(GPSummTable,2,as.character)
      write.table(
        GPSummTable,
        file = paste0(
          saveFilesTo,
          "\\ProActiveOutput\\ProActiveGenePredictstable.csv"
        ),
        sep = ",",
        row.names = FALSE
      )
    }
    write.table(
      summaryTable,
      file = paste0(
        saveFilesTo,
        "\\ProActiveOutput\\ProActiveSummaryTable.csv"
      ),
      sep = ",",
      row.names = FALSE
    )
    write.table(
      cleanSummaryTable,
      file = paste0(
        saveFilesTo,
        "\\ProActiveOutput\\ProActiveCleanSummaryTable.csv"
      ),
      sep = ",",
      row.names = FALSE
    )
    write.table(
      filteredOutContigsDf,
      file = paste0(
        saveFilesTo,
        "\\ProActiveOutput\\ProActiveFilteredOut.csv"
      ),
      sep = ",",
      row.names = FALSE
    )
    return(finalSummaryList)
  } else {
    return(finalSummaryList)
  }
}
