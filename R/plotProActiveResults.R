#' Plot results of `ProActive()` pattern-matching
#'
#' Plot read coverage of contigs/chunks with detected gaps and elevations and their
#' associated pattern-match.
#'
#' @param pileup  A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param qualityFilter Optional, only plot results with pattern-matches that achieved
#' a specified quality level. Options are "high", "okay", and "poor". Default is NA
#' (i.e. no filter).
#' @param ProActiveResults The output from `ProActive()`.
#' @param elevFilter Optional, only plot results with pattern-matches that achieved an
#' elevation ratio (max/min) greater than the specified values. Defualt is NA (i.e. no filter).
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#' output to. A folder will be made within the provided directory to store
#' results.
#' @import ggplot2
#' @return A list containing ggplot objects
#' @export
#' @examples
#' ProActivePlots <- plotProActiveResults(sampleMetagenomePileup,
#'                                        sampleMetagenomeResults)
plotProActiveResults <- function(pileup, ProActiveResults, qualityFilter = NA, elevFilter = NA, saveFilesTo) {
  position <- coverage <- NULL
  windowSize <- ProActiveResults[[5]][[1]]
  mode <- ProActiveResults[[5]][[2]]
  chunkSize <- ProActiveResults[[5]][[3]]
  contigChunk <- ProActiveResults[[5]][[4]]
  summaryTable <- ProActiveResults[[1]]
  patternMatches <- ProActiveResults[[3]]
  pileup <- pileupFormatter(pileup, mode)
  plots <- list()
  refNames <- rep(NA, nrow(summaryTable))
  if (mode == "genome") {
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if (contigChunk == "TRUE") {
    pileup <- contigChunks(pileup, chunkSize)
  }
  X <- 1
  lapply(seq_along(patternMatches), function(i) {
    refName <- patternMatches[[i]][[8]]
    matchInfo <- summaryTable[which(summaryTable[, 1] == refName), ]
    classification <- matchInfo[, 2]
    quality <- matchInfo[, 3]
    elevRatio <- ifelse(classification == "NoPattern", "NA", round(matchInfo[, 4], digits = 3))
    if (is.na(qualityFilter) == FALSE) {
      if (quality != qualityFilter) {
        return(NULL)
      }
    }
    if (is.na(elevFilter) == FALSE) {
      if (elevRatio < elevFilter) {
        return(NULL)
      }
    }
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    pattern <- patternBuilder(pileupSubset, patternMatches[[i]])
    patternMatch <- cbind(pileupSubset, pattern)
    matchLength <- matchInfo[, 7]
    plot <- ggplot(data = patternMatch, aes(x = position, y = coverage)) +
      geom_area(data = patternMatch, aes(x = position, y = coverage), fill = "seagreen3") +
      geom_line(y = pattern, linewidth = 1) +
      labs(title = paste(refName, classification), subtitle = paste("Matching-region size (bp):", matchLength, "Elevation ratio:", elevRatio), x = "contig position") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.subtitle = element_text(size = 12),
        plot.title = element_text(size = 14),
        plot.margin = margin(
          t = 0,
          r = 6,
          b = 0,
          l = 2
        )
      )
    plots[[X]] <<- plot
    refNames[[X]] <<- refName
    X <<- X + 1
  })
  refNames <- na.omit(refNames)
  names(plots) <- refNames
  if (missing(saveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\ProActiveOutput")),
      dir.create(paste0(saveFilesTo, "\\ProActiveOutput")),
      stop(
        "'ProActiveOutput' folder exists already in the provided
                directory"
      )
    )
    lapply(
      names(plots),
      function(X) {
        ggsave(
          filename = paste0(
            saveFilesTo,
            "\\ProActiveOutput\\", X, ".png"
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
