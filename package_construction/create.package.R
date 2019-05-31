# Commands to create MetaCycData package

library(devtools)
library(dplyr)

setwd("~/Files/MinT/github/MetaCycData")

# only once:
# setup(rstudio=FALSE)

# Load data from RDS files and save in Data directory

# data.dir <- "D:/Files/MinT/Metacyc biocyc files/meta/data"
data.dir <- "/Users/d3l348/tmp/meta-23.0/data"

mc_compounds <- readRDS(file.path(data.dir, "compounds.RDS"))
use_data(mc_compounds, overwrite=TRUE)

mc_reactions <- readRDS(file.path(data.dir, "reactions.RDS"))
use_data(mc_reactions, overwrite=TRUE)

mc_modules <- readRDS(file.path(data.dir, "modules.RDS"))
use_data(mc_modules, overwrite=TRUE)

# mc_pathways <- readRDS(file.path(data.dir, "pathways.RDS"))
# use_data(mc_pathways, overwrite=TRUE)


# Construct mappings between pathways, reactions and compounds

# reaction --> compound mapping
mc_reaction_compound_map <- strsplit(mc_reactions$COMPOUNDS, ";")
names(mc_reaction_compound_map) <- mc_reactions$REACTION
use_data(mc_reaction_compound_map, overwrite=TRUE)

# reaction --> module mapping
mc_module_reaction_map <- strsplit(mc_modules$`REACTION-LIST`, ";")
names(mc_module_reaction_map) <- mc_modules$MODULE
use_data(mc_module_reaction_map, overwrite=TRUE)

# compound --> reaction mapping
rxn.cmp.df <- do.call(rbind, lapply(1:length(mc_reaction_compound_map), function(i) {
    data.frame(reaction=names(mc_reaction_compound_map)[i], 
               compound=mc_reaction_compound_map[[i]], stringsAsFactors = FALSE)
}))
rxn.cmp.df <- group_by(rxn.cmp.df, compound) %>%
    tidyr::nest(reaction)
mc_compound_reaction_map <- rxn.cmp.df$data
names(mc_compound_reaction_map) <- rxn.cmp.df$compound
mc_compound_reaction_map <- lapply(mc_compound_reaction_map, function(x) as.vector(x$reaction))
    
use_data(mc_compound_reaction_map, overwrite=TRUE)

# for compound <--> formula mapping identify 1-many relationships
compound_formula_df <- select(mc_compounds, COMPOUND, MF)
mc_formulas_per_compound <- group_by(compound_formula_df, COMPOUND) %>%
  summarise(Num_Formula=n())
any(mc_formulas_per_compound$Num_Formula > 1) # FALSE

mc_compounds_per_formula <- group_by(compound_formula_df, MF) %>%
  summarise(Num_Compounds=n())
use_data(mc_compounds_per_formula, overwrite=TRUE)


# compound --> module node (1 or more reactions) mapping
cmp.rxn.module.df <- tibble::tibble(Compound=names(mc_compound_reaction_map), Reaction=mc_compound_reaction_map) %>%
  tidyr::unnest()
tmp <- tibble::tibble(Module=names(mc_module_reaction_map), Reaction=mc_module_reaction_map) %>%
  tidyr::unnest()
cmp.rxn.module.df <- inner_join(cmp.rxn.module.df, tmp) %>% unique()

# function to parse module nodes (1 or more reaction IDs) from mc_modules$REACTION field
getModuleNodes <- function(rxnString) {
  rxn.lines <- strsplit(rxnString, ";")[[1]]
  
  rxns <- character(0)
  nodes <- character(0)
  
  # remove lines that are just //
  ind.keep <- !grepl("^[//]+$", rxn.lines)
  rxn.lines <- rxn.lines[ind.keep]
  
  line.parts <- strsplit(rxn.lines, "(", fixed=TRUE)
  unlist(lapply(line.parts, function(x) x[2]))
  rxn.ids <- trimws(unlist(lapply(line.parts, function(x) x[2])))
  
  node.names <- data.frame(REACTION=rxn.ids, MODULE_NODE=rxn.ids, stringsAsFactors = FALSE)
  return(node.names)
}
tmp <- lapply(1:nrow(mc_modules), function(i) {
  .data <- getModuleNodes(mc_modules$`REACTION-LAYOUT`[i])
  .data$MODULE=mc_modules$MODULE[i]
  return(.data)
})
saveRDS(tmp, file=file.path(data.dir, "mc.rxn.module.node.mapping.RDS"))

con <- file(file.path(data.dir, "mc.rxn.module.node.mapping.csv"), "w")
cat(colnames(tmp[[1]]), sep=",", file=con)
cat("\n", file=con)
for (i in 1:length(tmp)) {
  write.table(tmp[[i]], col.names=FALSE, row.names=FALSE, sep=",", quote=FALSE, file=con)    
}
close(con)

mc_reaction_module_node_map <- read.table(file.path(data.dir, "mc.rxn.module.node.mapping.csv"), sep=",", header=TRUE, row.names=NULL)
mc_reaction_module_node_map <- unique(mc_reaction_module_node_map)

# make sure all reactions map to mc_reactions$REACTION
all(mc_reaction_module_node_map$REACTION %in% mc_reactions$REACTION)

use_data(mc_reaction_module_node_map, overwrite=TRUE)


## Get MetaCyc release number from version.dat file
version_info <- read.table(file.path(data.dir, "version.dat"), sep="\t", comment.char=";", row.names = 1)

message(sprintf("TODO: open DESCRIPTION file and increment version number and incorporate MetaCyc release info into title and description field: (Version %s, Released %s)", version_info["VERSION",1], version_info["RELEASE-DATE", 1]))


document()
# build()
install()
