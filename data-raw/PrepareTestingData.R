## Code to prepare datasets used for unit testing

library(ProActive)

## Genome
sampleGenomePileup <- read.delim("C:/Users/jlmaier/Desktop/ProActive/LT2_bincov100", header=FALSE, comment.char="#")
sampleGenomegffTSV <- read.delim("C:/Users/jlmaier/Desktop/ProActive/LT2.gff", header=FALSE, comment.char="#")

sampleGenomeResults <- ProActiveDetect(sampleGenomePileup, mode="genome", gffTSV = sampleGenomegffTSV)

usethis::use_data(sampleGenomeResults, internal = TRUE, overwrite = TRUE)

## Metagenome
M_spades <- read.delim("Q:/Shared drives/JessieMaier/Experiments/Transductomics/Data/Microbial_reads/M_spades.bincov100", header=FALSE, comment.char="#")

exp205WCcleaned <- read.delim("C:/Users/jlmaier/Desktop/ProActive/exp205WCcleaned.gff", header=FALSE)

pileupFormatter <- function(pileup) {
  colClasses <- vapply(pileup, class, character(1))
  for (i in c(which(colClasses == "integer"))) {
    if (length(which(pileup[, i] == 100)) > 1) {
      posColIdx <- i
    }
  }
  cleanPileup <-
    cbind.data.frame(
      pileup[, which(colClasses == "character")],
      pileup[, which(colClasses == "numeric")],
      pileup[, posColIdx]
    )
  colnames(cleanPileup) <- c("contigName", "coverage", "position")
  cleanPileup$contigName <- gsub("\\s.*", "", cleanPileup$contigName)
  return(cleanPileup)
}

M_spadesCleaned <- pileupFormatter(M_spades)

NODE_1911 <- M_spades[which(M_spadesCleaned[,1]=="NODE_1911"),]#NODE_1911: elevation off left
NODE_1583 <- M_spades[which(M_spadesCleaned[,1]=="NODE_1583"),]#NODE_1583: elevation off right
NODE_1884 <- M_spades[which(M_spadesCleaned[,1]=="NODE_1884"),]#NODE_1884: gap off right
NODE_1255 <- M_spades[which(M_spadesCleaned[,1]=="NODE_1255"),]#NODE_1255: gap off left
NODE_368 <- M_spades[which(M_spadesCleaned[,1]=="NODE_368"),]#NODE_368: full gap
NODE_617 <- M_spades[which(M_spadesCleaned[,1]=="NODE_617"),]#NODE_617: elevation full
NODE_1625 <- M_spades[which(M_spadesCleaned[,1]=="NODE_1625"),]#NODE_1625: no pattern

NODE_1911ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_1911"),]#NODE_1911: elevation off left
NODE_1583ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_1583"),]#NODE_1583: elevation off right
NODE_1884ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_1884"),]#NODE_1884: gap off right
NODE_1255ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_1255"),]#NODE_1255: gap off left
NODE_368ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_368"),]#NODE_368: full gap
NODE_617ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_617"),]#NODE_617: elevation full
NODE_1625ORFS <- exp205WCcleaned[which(exp205WCcleaned[,1]=="NODE_1625"),]#NODE_1625: no pattern

sampleMetagenomePileup <- rbind.data.frame(NODE_1911, NODE_1583, NODE_1884,
                                           NODE_1255, NODE_368, NODE_617, NODE_1625)

sampleMetagenomegffTSV <- rbind.data.frame(NODE_1911ORFS, NODE_1583ORFS, NODE_1884ORFS,
                                           NODE_1255ORFS, NODE_368ORFS, NODE_617ORFS, NODE_1625ORFS)

sampleMetagenomeResults <- ProActiveDetect(sampleMetagenomePileup, mode="metagenome", gffTSV = sampleMetagenomegffTSV)

usethis::use_data(sampleMetagenomeResults, overwrite = TRUE)
