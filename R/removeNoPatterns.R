
removeNoPatterns <- function(bestMatchList){
  newBestMatchList <- list()
  length(bestMatchList)
  X<-1
  lapply(seq_along(bestMatchList), function(i){
    bestMatchInfo <- bestMatchList[[i]]
    classification <- bestMatchInfo[[7]]
    if (classification == "NoPattern") {
      return (NULL)
    } else {
      newBestMatchList[[X]] <<- bestMatchInfo
      X<<-X+1
    }
  })

  return(newBestMatchList)
}
