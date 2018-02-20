#' Mapping from reactions to compounds
#' 
#' A list that contains a vector of compound names for each reaction. 
#' 
#' @format A named list of length 15311, where the names correspond to 
#' \code{mc_reaction$REACTION} and the values are a list of compound names 
#' corresponding to \code{mc_compounds$COMPOUND}.
#' @source \url{http://metacyc.ai.sri.com/}
#' @seealso \code{\link{mc_reactions}}, \code{\link{mc_compounds}}
"mc_reaction_compound_map"