#' Mapping from modules to reactions
#' 
#' A list that contains a vector of reaction IDs for each module. 
#' 
#' @format A named list of length 2572, where the names correspond to 
#' \code{mc_modules$MODULE} and the values are a list of reaction IDs 
#' corresponding to \code{mc_reactions$REACTION}.
#' @source \url{http://metacyc.ai.sri.com/}
#' @seealso \code{\link{mc_reactions}}, \code{\link{mc_modules}}
"mc_module_reaction_map"