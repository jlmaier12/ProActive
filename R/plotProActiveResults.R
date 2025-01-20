#' Plot results of `ProActive()` pattern-matching
#'
#' Plot read coverage of contigs/chunks with detected gaps and elevations and their
#' associated pattern-match.
#'
#' @param pileup  A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param ProActiveResults The output from `ProActive()`.
#' @param elevFilter Optional, only plot results with pattern-matches that achieved an
#' elevation ratio (max/min) greater than the specified values. Default is no filter.
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#' output to. A folder will be made within the provided directory to store
#' results.
#' @import ggplot2
#' @return A list containing ggplot objects
#' @export
#' @examples
#' ProActivePlots <- plotProActiveResults(sampleMetagenomePileup,
#'                                        sampleMetagenomeResults)
plotProActiveResults <- function(pileup, ProActiveResults, elevFilter, saveFilesTo) {
  position <- coverage <- pattern <- NULL
  windowSize <- ProActiveResults[[5]][[1]]
  mode <- ProActiveResults[[5]][[2]]
  chunkSize <- ProActiveResults[[5]][[3]]
  contigChunk <- ProActiveResults[[5]][[4]]
  summaryTable <- ProActiveResults[[1]]
  patternMatches <- ProActiveResults[[3]]
  pileup <- pileupFormatter(pileup, mode)
  if (mode == "genome") {
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if (mode == "metagenome" & contigChunk == TRUE) {
    pileup <- contigChunks(pileup, chunkSize)
  }
  elevFilter <- ifelse(missing(elevFilter), NA, elevFilter)
  plots <- lapply(seq_along(patternMatches), function(i) {
    refName <- patternMatches[[i]][[8]]
    matchInfo <- summaryTable[which(summaryTable[, 1] == refName), ]
    classification <- matchInfo[, 2]
    elevRatio <- ifelse(classification == "NoPattern", 0, round(matchInfo[, 3], digits = 3))
    if (is.na(elevFilter) == FALSE & elevRatio < elevFilter) {
      return(NULL)
    }
    pileupSubset <- pileup[which(pileup[, 1] == refName), ]
    pileupSubset <- changewindowSize(pileupSubset, windowSize, mode)
    patternMatch <- patternBuilder(pileupSubset, patternMatches[[i]])
    matchLength <- matchInfo[, 6]
    plot <- ggplot(data = patternMatch, aes(x = position, y = coverage)) +
      geom_area(fill = "seagreen3") +
      geom_line(data = patternMatch, aes(y = pattern), linewidth = 1) +
      labs(title = paste(refName, classification),
           subtitle = paste("Matching-region size (bp):", matchLength, "Elevation ratio:", elevRatio),
           x = "Base pair (bp) position") +
      scale_x_continuous(expand = c(0, 0)) +
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
          r = 10,
          b = 0,
          l = 2
        )
      )
    plot
  })
  plots <- (plots[!vapply(plots, is.null, logical(1))])
  refNames <- sapply(seq_along(patternMatches), function(i) {
    refName <- patternMatches[[i]][[8]]
    matchInfo <- summaryTable[which(summaryTable[, 1] == refName), ]
    classification <- matchInfo[, 2]
    elevRatio <- ifelse(classification == "NoPattern", 0, round(matchInfo[, 3], digits = 3))
    if (is.na(elevFilter) == FALSE & elevRatio < elevFilter) {return(NA)}
    refName
  })
  refNames <- (refNames[!vapply(refNames, is.na, logical(1))])
  names(plots) <- refNames
  if (missing(saveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\ProActiveOutputPlots")),
      dir.create(paste0(saveFilesTo, "\\ProActiveOutputPlots")),
      stop(
        "'ProActiveOutputPlots' already exists in the provided directory"
      )
    )
    lapply(
      names(plots),
      function(X) {
        ggsave(
          filename = paste0(
            saveFilesTo,
            "\\ProActiveOutputPlots\\", X, ".png"
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
