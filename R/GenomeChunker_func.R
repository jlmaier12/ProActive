#' Break genome file into chunks for pattern-matching analysis
#'
#' @param genome_pileup Pileup file for read coverages mapped to a genome
#' @param chunk_size Value to chunk genome into. Should be divisible by 100. Default is 100,000bp per chunk.
#'
#' @keywords internal
GenomeChunks <- function(genome_pileup, chunk_size) {
  ref_name <- rep(NA, nrow(genome_pileup))
  position <- rep(NA, nrow(genome_pileup))
  genome_ref <- genome_pileup[1,1]
  chunk_size <- chunk_size/100
  X <- 1
  Y <- chunk_size
  Z <- 1
  repeat {
  ref_name[c(X:Y)] <- rep(paste0(genome_ref,"_chunk_", Z), chunk_size)
  if(Z==max(genome_pileup$position)%/%(chunk_size*100)) break
  X <- X+chunk_size
  Y <- Y+chunk_size
  Z <- Z+1
  }
  Z <- Z+1
  remain_index <- which(is.na(ref_name))[1]
  genome_chunk <- genome_pileup[remain_index:nrow(genome_pileup),]
  ref_name[c(remain_index:nrow(genome_pileup))] <- rep(paste0("genome_chunk", Z), nrow(genome_chunk))
  chunked_genome <- cbind.data.frame(ref_name, genome_pileup$coverage, genome_pileup$position)
  colnames(chunked_genome) <- c("ref_name", "coverage", "position")
  return(chunked_genome)
}


