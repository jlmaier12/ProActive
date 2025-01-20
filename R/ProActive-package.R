#' @title \bold{ProActive}
#'
#' @description `ProActive` automatically detects regions of gapped and elevated read coverage
#' using a 2D pattern-matching algorithm. `ProActive` detects, characterizes and
#' visualizes read coverage patterns in both genomes and metagenomes. Optionally,
#' users may provide gene annotations associated with their genome or metagenome
#' in the form of a .gff file. In this case, `ProActive` will generate an additional
#' output table containing the gene annotations found within the detected regions of
#' gapped and elevated read coverage. Additionally, users can search for gene
#' annotations of interest in the output read coverage plots.
#'
#' @details The three main functions in `ProActive` are:
#' \enumerate{
#' \item \code{\link{ProActiveDetect}} performs the pattern-matching and characterization of read coverage patterns.
#' \item \code{\link{plotProActiveResults}} plots the results from
#'  \code{ProActiveDetect()}
#' \item \code{\link{geneAnnotationSearch}} searches classified contigs/chunks for gene annotations that match
#' user-provided keywords.
#' }
#'
#' @author Jessie Maier \email{jlmaier@ncsu.edu}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
