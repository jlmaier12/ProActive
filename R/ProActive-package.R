#' @title \bold{ProActive}
#'
#' @description Automatic detection of gaps and elevations in mapped sequencing reads using read coverage
#'   pattern-matching. Pattern-matching results are summarized into a dataframe and users can visualize the
#'   pattern-matches and associated read coverage patterns. If a .gff file is supplied, ProActive will output
#'   the ORFs found within the detected regions of elevated and/or gapped read coverage. ProActive works with both
#'   metagenomes and genomes.
#'
#' @details The two main functions in ProActive are:
#' \enumerate{
#' \item \code{\link{ProActive}} performs the pattern-matching and characterization of read coverage patterns.
#' \item \code{\link{plotProActiveResults}} plots the results from
#'  \code{ProActive()}
#' }
#'
#' @author Jessie Maier \email{jlmaier@ncsu.edu}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
