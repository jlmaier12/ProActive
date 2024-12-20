#' Full elevation pattern pattern-translator for whole-community read coverage only
#'
#' Translates full elevation pattern across a contig. Translate the pattern 1000bp at a time. Stop translating when the elevation pattern is 5000bp from the end of the contig.
#'
#' @param contigCov The read coverages that pertains only to the contig currently being assessed
#' @param bestMatchInfo The information associated with the current best pattern match. Includes the comparison score, the minimum and maximum pattern y-acis values, the pattern length, and the pattern start and stop positions on the contig
#' @param windowSize The window size used to re-average read coverage datasets
#' @param pattern A vector containing the values associated with the elevation pattern to be translated across the contig
#' @param elevOrGap
#'
#' @keywords internal
patternTranslator <- function(contigCov, bestMatchInfo, windowSize, pattern, elevOrGap){
  minPattern <- min(pattern)
  maxPattern <- max(pattern)
  minOrMax <- ifelse(elevOrGap=="Elevation", min(pattern), max(pattern))
  maxOrMin2 <- ifelse(elevOrGap=="Elevation", max(pattern), min(pattern))
  repeat {
    pattern <- c(rep(minOrMax, (1000/windowSize)),pattern[-c((length(pattern)-((1000/windowSize)-1)):length(pattern))])
    if(elevOrGap=="Elevation"){
      if(pattern[length(pattern)-(5000/windowSize)] > minPattern) break
      } else {
      if(pattern[length(pattern)-(5000/windowSize)] < maxPattern) break
    }
    diff <- mean(abs(contigCov - pattern))
    if (diff < bestMatchInfo[[1]]){
      elevRatio <- max(pattern)/min(pattern)
      startPos <- which(pattern==maxOrMin2)[1]
      endPos <- which(pattern==maxOrMin2)[length(which(pattern==maxOrMin2))]
      bestMatchInfo <- list(diff, minPattern, maxPattern, startPos, endPos, elevRatio, elevOrGap)
    }
  }
  return(bestMatchInfo)
}
