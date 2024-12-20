#' Plot read coverage graphs of contigs with predicted active prophages in the whole-community read coverages
#'
#' Plot the read coverages of a contig and its associated pattern-match for each contig with a predicted active prophage.
#'
#' @param pileup A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param qualityFilter A list containing pattern information associated with all contigs containing a potential active prophage. The second, third or fourth item in the list output from ProActive()
#' @param ProActiveResults The result from `ProActive()`
#' @param elevFilter Default is NA. The user may define an elevation factor in which any plots with elevated read coverages below the elevation factor are not displayed
#' @param saveFilesTo
#' @export
#' @examples
#' \dontrun{
#' data(samplePileup)
#'
#' plotProActiveResults(whole_commreadcovs, ProActiveResults, ProActiveResults$ExtraConfidentclassifications)
#' }
plotProActiveResults <- function(pileup, ProActiveResults, qualityFilter=NA, elevFilter=NA, saveFilesTo) {
  position <- coverage <- NULL
  windowSize <- ProActiveResults[[4]][[1]]
  mode <- ProActiveResults[[4]][[2]]
  chunkSize <- ProActiveResults[[4]][[3]]
  contigChunk <- ProActiveResults[[4]][[4]]
  summaryTable <- ProActiveResults[[1]]
  patternMatches <- ProActiveResults[[2]]
  pileup <- pileupFormatter(pileup, mode)
  plots <- list()
  refNames <- rep(NA, nrow(summaryTable))
  if (mode=="genome"){
    pileup <- genomeChunks(pileup, chunkSize)
  }
  if(contigChunk =="TRUE"){
    pileup <- contigChunks(pileup, chunkSize)
  }
  X <- 1
  lapply(seq_along(patternMatches),function(i){
    refName <- patternMatches[[i]][[8]]
    matchInfo <- summaryTable[which(summaryTable[,1]==refName),]
    classification <- matchInfo[,2]
    quality <- matchInfo[,3]
    elevRatio <- ifelse(classification=="NoPattern","NA",round(matchInfo[,4], digits=3))
    if(is.na(qualityFilter)==FALSE){
      if(quality != qualityFilter){
        return(NULL)
      }
    }
    if(is.na(elevFilter)==FALSE){
      if(elevRatio < elevFilter){
        return(NULL)
      }
    }
    pileupSubset <- pileup[which(pileup[,1] == refName),]
    pileupSubset <- changewindowSize(pileupSubset,windowSize, mode)
    pattern <- patternBuilder(pileupSubset, patternMatches[[i]])
    patternMatch <- cbind(pileupSubset, pattern)
    matchLength <- matchInfo[,7]
    plot <- ggplot(data=patternMatch, aes(x=position, y=coverage))+
            geom_area(data=patternMatch, aes(x=position, y=coverage), fill="seagreen3") +
            geom_line(y=pattern, linewidth=1)+
            labs(title=paste(refName, classification), subtitle=paste("Matching-region size (bp):", matchLength, "Elevation ratio:", elevRatio), x="contig position") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  text = element_text(size = 15))
    plots[[X]] <<- plot
    refNames[[X]] <<- refName
    X<<-X +1
  })
  refNames <- na.omit(refNames)
  names(plots) <- refNames
  if(missing(saveFilesTo)==FALSE){
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
