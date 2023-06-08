#' Prokka gff file reformatter
#'
#' Prokka gene annotations cleanup
#'
#' @param gff_file A gff file containing gene annotation information from Prokka. Your gff file should have columns in the following order: contig (i.e "NODE_#"), tool (i.e "prodigal"), CDS, start position, stop position, '.', '+/-', '0', annotations.
#'
#' @export
#'
prokka_cleanup <- function(gff_file) {
  V9 <- NULL
  gff_df <- gff_file[-c(2,3,6,7,8)] #***need to fix this so correct columns are choosen every time
  colnames(gff_df)[1] <- "ref_name"
  gff_df <- separate(gff_df, col=V9, into=c('ID', 'annotations'), sep=';', extra="merge")
  inference <- str_extract(gff_df[,5], "(?<=inference=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  EC_number <- str_extract(gff_df[,5], "(?<=eC_number=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  gene <- str_extract(gff_df[,5], "(?<=gene=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  product <- str_extract(gff_df[,5], "(?<=;product=)[\\s\\S]*")
  cleaned_gff_df <- cbind(gff_df[,c(1:4)], inference, EC_number, gene, product)
  cleaned_gff_df[,1] <- gsub("_length.*", "", gff_df[,1])
  colnames(cleaned_gff_df)[c(1:4)] <- c("ref_name", "startpos", "stoppos","annotation")
  return(cleaned_gff_df)
}
