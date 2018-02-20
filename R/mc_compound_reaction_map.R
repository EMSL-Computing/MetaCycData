#' Mapping from compounds to reactions
#' 
#' A list that contains a vector of reaction names for each compound. 
#' 
#' @format A named list of length 12875, where the names correspond to 
#' \code{mc_compounds$COMPOUND} and the values are a list of compound names 
#' corresponding to \code{mc_reaction$REACTION}.
#' @source \url{http://metacyc.ai.sri.com/}
#' @seealso \code{\link{mc_reactions}}, \code{\link{mc_compounds}}
"mc_compound_reaction_map"