# Parse MetaCyc pathways.dat file
library(dplyr)

data.dir <- "/Users/d3l348/tmp/meta-23.0/data"
fname <- file.path(data.dir, "pathways.dat")

# Find the starts of fields (i.e. all caps title at beginning of line)
field.start <- function(txt) {
    ind <- grep("^[[:upper:]]+", txt)
    return(ind)
}

field.type <- function(txt) {
    ind <- regexpr("^[[:upper:]_-]+", txt)
    fields <- unlist(lapply(1:length(ind), 
                            function(i) { ifelse(ind[i] > 0, 
                                                 substring(txt[i], ind[i], ind[i]+attr(ind, "match.length")[i]-1), NA) 
                            }))
    fields <- fields[!is.na(fields)]
    return(fields)
}

# Returns the field value including any embedded returns
# 'txt' should be only one field text, not multiple
field.value <- function(txt) {
    ftype <- field.type(txt)
    #remove field name
    txt <- gsub(paste0("^",ftype), "", txt)
    #remove ' - ' after field name
    txt <- gsub("^ - ", "", txt)
    #de-indent
    txt <- gsub("^[[:space:]]*", "", txt)
    #remove trailing spaces
    txt <- gsub("[[:space:]]*$", "", txt)
    #remove end of compound marker /// and empty lines
    txt <- gsub("//", "", txt)
    txt <- txt[nchar(txt) > 0]
    #collapse multiple lines
    txt <- paste(txt, collapse="\n")
    txt <- gsub(" *\n", "\n", txt)
    return(txt)
}


database <- scan(fname, what = "character", sep="\n", fileEncoding = "windows-1252")

pth.sts <- grep("^UNIQUE-ID", database)
pth.ends <- c(pth.sts[-1] - 1, length(database))

res <- lapply(1:length(pth.ends), function(x) NULL)

for(cur in 1:length(pth.ends)){ 
    
    cur.pth <- database[pth.sts[cur]:pth.ends[cur]]
    
    field.sts <- field.start(cur.pth)
    field.ends <- c(field.sts[-1]-1, length(cur.pth))
    
    field.names <- c("UNIQUE-ID", "URL", "COMMON-NAME", "TYPES", "SUPER-PATHWAYS",  "SUB-PATHWAYS", 
                    "REACTION-LIST", "REACTION-LAYOUT", "SPECIES", "TAXONOMIC-RANGE", "PREDECESSORS", 
                    "PATHWAY-LINKS", "COMMENT", "DBLINKS")
    
    # Currently ignoring fields: EXACT, MOL, MODULE
    field.values <- rep(NA, times=length(field.names))
    names(field.values) <- field.names
    
    for (ff in 1:length(field.sts)) {
        cur.field <- cur.pth[field.sts[ff]:field.ends[ff]]
        ftype <- field.type(cur.field)
        fvalue <- field.value(cur.field)
        
        #any other adaptation of field text needed before saving?
        if (ftype == "UNIQUE-ID") {
            fvalue <- strsplit(fvalue, " +")[[1]][1]
            mc.url <- paste("http://metacyc.ai.sri.com/META/new-image?type=PATHWAY&object=",fvalue, sep = "")
            field.values["URL"] <- mc.url
        } else if (ftype == "COMMENT" & grepl("\n", fvalue)) {
            #browser()
            fvalue <- gsub("\n/", "\n", fvalue)
        } 
        
        if (ftype %in% field.names & is.na(field.values[ftype])) {
            field.values[ftype] <- fvalue
        } else if (ftype %in% c("REACTION-LIST", "REACTION-LAYOUT", "PREDECESSORS")) {
            field.values[ftype] <- paste(fvalue, field.values[ftype], sep=";")
        } else if (ftype %in% field.names) {
            field.values[ftype] <- paste(field.values[ftype], fvalue, sep=";")
        }
        
    }
    
    res[[cur]] <- field.values
}	

# Turn into a data.frame
pathway_res <- do.call(rbind, res)

pathway_res <- as.data.frame(pathway_res, stringsAsFactors=FALSE)

field.names <- gsub("UNIQUE-ID", "PATHWAY", field.names)
names(pathway_res) = field.names

pathway_res$IS_SUPER_PATHWAY <- grepl("Super-Pathways", pathway_res$TYPES)
modules <- filter(pathway_res, !IS_SUPER_PATHWAY)
modules <- rename(modules, MODULE=PATHWAY)

pathway_res <- filter(pathway_res, IS_SUPER_PATHWAY)
saveRDS(pathway_res, file.path(data.dir, "pathways.RDS"))

saveRDS(modules, file.path(data.dir, "modules.RDS"))
