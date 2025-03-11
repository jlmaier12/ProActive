#' Gene annotation plot
#'
#' Plot read coverage and location of gene annotations that match the keywords and
#' search criteria for contig/chunk currently being assessed
#'
#' @param geneAnnotSubset Subset of gene annotations to be plotted
#' @param keywords The key-word(s) used for the search.
#' @param pileupSubset A subset of the pileup associated with the contig/chunk being assessed
#' @param colIdx The column index 'gene' or 'product' column
#' @param startbpRange The basepair at which the search is started if a 'specific' search is used
#' @param endbpRange The basepair at which the search is ended if a 'specific' search is used
#' @param elevRatio The maximum/minimum values of the pattern-match
#' @param pattern The pattern-match information associated with the contig/chunk being assessed
#' @param windowSize The number of basepairs to average read coverage values over.
#' @param chunkSize If `mode`="genome" OR if `mode`="metagenome" and `chunkContigs`=TRUE,
#' chunk the genome or contigs, respectively, into smaller subsets for pattern-matching.
#' `chunkSize` determines the size (in bp) of each 'chunk'. Default is 50000.
#' @param mode Either "genome" or "metagenome"
#' @keywords internal
#' @importFrom stringr str_which
geneAnnotationPlot <- function(geneAnnotSubset, keywords, pileupSubset,
                               colIdx, startbpRange, endbpRange, elevRatio,
                               pattern, windowSize, chunkSize, mode, contigChunk){
  position <- coverage <- start <- NULL
  classification <- pattern[[7]]
  refName <- pileupSubset[1, 1]
  if(mode == "genome"){
    chunkNumber <- as.numeric(str_extract(refName, "(?<=\\_)\\d+$")) - 1
    startPos <- (pattern[[4]] * windowSize) + (chunkNumber * chunkSize)
    endPos <- (pattern[[5]] * windowSize) + (chunkNumber * chunkSize)
  } else if (grepl("chunk", refName, fixed = TRUE)){
    chunkNumber <- as.numeric(str_extract(refName, "(?<=\\_)\\d+$")) - 1
    startPos <- (pattern[[4]] * windowSize) + (chunkNumber * chunkSize)
    endPos <- (pattern[[5]] * windowSize) + (chunkNumber * chunkSize)
  }else {
    startPos <- pattern[[4]] * windowSize
    endPos <- pattern[[5]] * windowSize
  }
  matchIdxs <- str_which(geneAnnotSubset[, colIdx], regex(paste(keywords, collapse = "|"), ignore_case = TRUE))
  geneAnnotMatches <- geneAnnotSubset[matchIdxs,]
  geneStartPos <- geneAnnotMatches$start
  geneAnnotLabels <- paste0("#", c(1:nrow(geneAnnotMatches)), ": ", geneAnnotMatches[, colIdx],
                            sep = " ", collapse = " \n ")
  plot <- ggplot(data = pileupSubset, aes(x = position, y = coverage)) +
    geom_area(fill = "#009E73") +
    geom_vline(xintercept = geneStartPos, linewidth = 1) +
    geom_vline(xintercept = c(startPos, endPos), color = "#D55E00", linewidth = 1) +
    geom_vline(xintercept = c(startbpRange, endbpRange), color = "#D55E00", linewidth = 1, linetype = "dotted") +
    geom_label(data = geneAnnotMatches, aes(x = start, y = (max(pileupSubset$coverage) / 2),
                                           label = paste0("#", c(1 : nrow(geneAnnotMatches)))),
               size = 2.75) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 15),
          plot.margin = margin(
            t = 0,
            r = 10,
            b = 0,
            l = 2
          )) +
    labs(title = paste(refName, classification),
         subtitle = paste("elevation ratio:", round(elevRatio, digits = 4)),
         x = "Basepair position",
         caption = geneAnnotLabels,
         y = "Read coverage")
  return(plot)
}
