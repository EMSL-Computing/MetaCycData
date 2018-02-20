#' Module information from MetaCyc
#' 
#' A data.frame containing module information parsed from the MetaCyc database.
#' Note that in MetaCyc, these are called Pathways, however, "Modules" only include
#' ones that are not super-pathways.
#' 
#' @format A data frame with 2572 rows and 14 variables:
#' \describe{
#'   \item{MODULE}{ }
#'   \item{URL}{ }
#'   \item{COMMON-NAME}{ }
#'   \item{TYPES}{ }
#'   \item{SUPER-PATHWAYS}{ }
#'   \item{SUB-PATHWAYS}{ }
#'   \item{REACTION-LIST}{ }
#'   \item{SPECIES}{ }
#'   \item{TAXONOMIC-RANGE}{ }
#'   \item{PREDECESSORS}{ }
#'   \item{PATHWAY-LINKS}{ }
#'   \item{COMMENT}{ }
#'   \item{DBLINKS}{ }
#'   \item{IS_SUPER_PATHWAY}{ }
#' }
#' @source \url{http://metacyc.ai.sri.com/}
"mc_modules"