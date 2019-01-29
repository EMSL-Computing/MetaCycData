# Parse MetaCyc reactions.dat file


data.dir <- "/Users/d3l348/tmp/MetaCyc_local_copy/22.6/data"
fname <- file.path(data.dir, "reactions.dat")

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

rxn.sts <- grep("^UNIQUE-ID", database)
rxn.ends <- c(rxn.sts[-1] - 1, length(database))

res <- lapply(1:length(rxn.ends), function(x) NULL)

for(cur in 1:length(rxn.ends)){ 
    
    cur.comp <- database[rxn.sts[cur]:rxn.ends[cur]]
    
    field.sts <- field.start(cur.comp)
    field.ends <- c(field.sts[-1]-1, length(cur.comp))
    
    field.names <- c("UNIQUE-ID", "URL", "SYSTEMATIC-NAME", "EC-NUMBER", "IN-PATHWAY", 
                     "REACTION-LIST", "LEFT", "RIGHT", "REACTION-DIRECTION", 
                     "COMPOUNDS", "DBLINKS")
    
    # Currently ignoring fields: EXACT, MOL, MODULE
    field.values <- rep(NA, times=length(field.names))
    names(field.values) <- field.names
    
    for (ff in 1:length(field.sts)) {
        cur.field <- cur.comp[field.sts[ff]:field.ends[ff]]
        ftype <- field.type(cur.field)
        fvalue <- field.value(cur.field)
        
        #any other adaptation of field text needed before saving?
        if (ftype == "UNIQUE-ID") {
            fvalue <- strsplit(fvalue, " +")[[1]][1]
            mc.url <- paste("http://metacyc.ai.sri.com/META/new-image?type=REACTION&object=",fvalue, sep = "")
            field.values["URL"] <- mc.url
        } else if ((ftype == "LEFT" | ftype == "RIGHT") & grepl("COEFFICIENT", fvalue)) {
            lines <- unlist(strsplit(fvalue, "\n"))
            coef <- gsub("\\^COEFFICIENT *- *", "", lines[2])
            fvalue <- paste0("(", coef, " ", lines[1], ")")
            field.values["COMPOUNDS"] <- ifelse(is.na(field.values["COMPOUNDS"]), lines[1], 
                                                     paste(field.values["COMPOUNDS"], lines[1], sep=";"))
        } else if ((ftype == "LEFT" | ftype == "RIGHT") & grepl("COMPARTMENT", fvalue)) {
            lines <- unlist(strsplit(fvalue, "\n"))
            comp <- gsub("\\^COMPARTMENT *- *", "", lines[2])
            fvalue <- paste0(lines[1], " (COMPARTMENT ", comp, ")")
            field.values["COMPOUNDS"] <- ifelse(is.na(field.values["COMPOUNDS"]), lines[1], 
                                                paste(field.values["COMPOUNDS"], lines[1], sep=";"))
        } else if (ftype %in% c("RIGHT", "LEFT")) {
            field.values["COMPOUNDS"] <- ifelse(is.na(field.values["COMPOUNDS"]), fvalue, 
                                                     paste(field.values["COMPOUNDS"], fvalue, sep=";"))
        } else if (ftype == "EC-NUMBER" & grepl("\n", fvalue)) {
            fvalue <- unlist(strsplit(fvalue, "\n"))[1]
            fvalue <- gsub("^EC-", "", fvalue)
        } else if (ftype == "EC-NUMBER") {
            fvalue <- gsub("^EC-", "", fvalue)
        } 
        
        if (ftype %in% field.names & is.na(field.values[ftype])) {
            field.values[ftype] <- fvalue
        } else if (ftype %in% field.names & ftype %in% c("LEFT", "RIGHT")) {
            field.values[ftype] <- paste(fvalue, "+", field.values[ftype])
        } else if (ftype %in% field.names & ftype %in% c("IN-PATHWAY", "REACTION-LIST")) {
            field.values[ftype] <- paste(field.values[ftype], fvalue, sep=";")
        } else if (ftype %in% field.names) {
            field.values[ftype] <- paste(field.values[ftype], fvalue)
        }
        
    }
    
    res[[cur]] <- field.values
}	

# Turn into a data.frame
reaction_res <- do.call(rbind, res)

reaction_res <- as.data.frame(reaction_res, stringsAsFactors=FALSE)

field.names <- gsub("UNIQUE-ID", "REACTION", field.names)
names(reaction_res) = field.names

saveRDS(reaction_res, file.path(data.dir, "reactions.RDS"))
