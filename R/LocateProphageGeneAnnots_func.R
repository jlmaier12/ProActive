#' Locate phage-related gene annotations on contigs containing potential active prophages
#'
#' Search contigs containing potential active prophages (detected in whole-community read coverage only) for gene-annotations containing the words: "phage","bacteriophage","portal","spike", "tail", "capsid", "holin", "integrase", "needle", and "lysin". The matching gene-annotation locations are plotted on the associated read coverage graph
#'
#' @param cleaned_gene_annots A table of gene annotations made by PROKKA that have been passed through the prokka_cleanup function.
#' @param prophagepredictions A list containing shape information associated with all contigs containing a potential active prophage. The third item in the list output from the main prophage_activity_finder_func
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' phage_gene_finder(gene_annots_cleaned, active_prophage_results[[3]], whole_commreadcoverages)
#' }
phage_gene_finder <- function(cleaned_gene_annots,prophagepredictions, microbialread_dataset) {
  position <- coverage <- startpos <- NULL
  microbialread_dataset <- readcovdf_formatter(microbialread_dataset)
  for (i in  seq(1, length(prophagepredictions), 1)) {
    ref_name <- prophagepredictions[[i]][[7]]
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == ref_name),]
    gene_annots <- cleaned_gene_annots[which(cleaned_gene_annots[,1]== ref_name),]
    start_pos <- prophagepredictions[[i]][[4]]
    end_pos <- prophagepredictions[[i]][[5]]
    phage_genes <- c("phage","bacteriophage","portal","spike", "tail", "capsid", "holin", "integrase", "needle",  "lysin")
    if (TRUE%in% (str_detect(gene_annots$product, regex(paste(phage_genes, collapse="|"),ignore_case=T)))==TRUE) {
      match_indexes <- str_which(gene_annots$product, regex(paste(phage_genes, collapse="|"), ignore_case=T))
      if ((TRUE%in%(match_indexes %in% seq(start_pos, end_pos, 1)))==TRUE) {
        gene_annot_subset <- gene_annots[match_indexes,]
        gene_start_pos <- gene_annot_subset[,2]
        gene_annots_labels <- paste0("#", c(1:nrow(gene_annot_subset)), ": ",gene_annot_subset$product, sep=" ", collapse=" \n ")

        print(ggplot(data=microbial_subset, aes(x=position, y=coverage))+
                geom_area(fill="deepskyblue3") +
                geom_vline(xintercept=gene_start_pos, size=1)+
                geom_label(data=gene_annot_subset, aes(x=startpos,y=(max(microbial_subset$coverage)/2),label=paste0("#", c(1:nrow(gene_annot_subset)))), size=2.5)+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))+
                labs(title=ref_name, x=NULL, caption=gene_annots_labels))
      }
    }
  }
}
