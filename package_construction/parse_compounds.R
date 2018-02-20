# Parse MetaCyc compounds.dat file


# data.dir <- "D:/Files/MinT/Metacyc biocyc files/meta/data"
data.dir <- "~/Files/MinT/Data/MetaCyc_21.1/data"
fname <- file.path(data.dir, "compounds.dat")

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

comp.sts <- grep("^UNIQUE-ID", database)
comp.ends <- c(comp.sts[-1] - 1, length(database))

res <- lapply(1:length(comp.ends), function(x) NULL)

for(cur in 1:length(comp.ends)){ 
    
    cur.comp <- database[comp.sts[cur]:comp.ends[cur]]
    
    field.sts <- field.start(cur.comp)
    field.ends <- c(field.sts[-1]-1, length(cur.comp))
    
#     field.names <- c("ENTRY", "URL", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT", 
#                      "REACTION", "ENZYME", "PATHWAY", "COMMENT", "REMARK", "BRITE", "ATOM", "BOND", "DBLINKS")
    field.names <- c("UNIQUE-ID", "URL", "COMMON-NAME", "MONOISOTOPIC-MW", "MOLECULAR-WEIGHT", 
                     "MF", "CHEMICAL-FORMULA", "DBLINKS")
    
    
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
            kegg.url <- paste("http://metacyc.ai.sri.com/compound?orgid=META&id=",fvalue, sep = "")
            field.values["URL"] <- kegg.url
        }
        
        if (ftype %in% field.names & is.na(field.values[ftype])) {
            field.values[ftype] <- fvalue
        } else if (ftype %in% field.names & !is.na(field.values[ftype])) {
            field.values[ftype] <- paste(field.values[ftype], fvalue)
        }
        
    }
    
    res[[cur]] <- field.values
}	

# Turn into a data.frame
compound_res <- do.call(rbind, res)

compound_res <- as.data.frame(compound_res, stringsAsFactors=FALSE)

field.names <- gsub("UNIQUE-ID", "COMPOUND", field.names)
names(compound_res) = field.names

compound_res$`MONOISOTOPIC-MW` <- as.numeric(compound_res$`MONOISOTOPIC-MW`)
compound_res$`MOLECULAR-WEIGHT` <- as.numeric(compound_res$`MOLECULAR-WEIGHT`)

# calculate formula from CHEMICAL-FORMULA field
## Order: C H O N S P
df.cnames <- c("C", "H", "O", "N", "S", "P")

# first get all unique letters used
unique.letters <- unique(unlist(strsplit(gsub("[^[:alpha:]]+", " ", compound_res$`CHEMICAL-FORMULA`), " ")))
unique.letters <- sort(setdiff(unique.letters, c(df.cnames, "", NA)))
# convert two letter abbrevs to upper-lower form
.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
        sep = "", collapse = " ")
}
# unique.letters <- unlist(lapply(unique.letters, .simpleCap))
df.cnames <- c(df.cnames, unique.letters)
parseChemicalFormulaField <- function(cfVec) {
    parts <- strsplit(cfVec, "\\) \\(")
    molList <- lapply(parts, function(x) {
        counts <- rep(0, times=length(df.cnames) )
        names(counts) <- df.cnames
        
        x <- gsub("\\)$", "", gsub("^\\(", "", x))
        res <- matrix(unlist(lapply(x, function(y) {
            unlist(strsplit(y, " "))
        })), byrow = TRUE, ncol=2)
        
        counts <- counts+as.numeric(res[match(df.cnames, res[,1]), 2])
        counts[is.na(counts)] <- 0
        return(counts)
    })
    
    molMat <- do.call(rbind, molList)
    colnames(molMat) <- unlist(lapply(colnames(molMat), .simpleCap))
    mf <- apply(molMat, 1, FUN=function(x) {
      x <- x[x>0]
      return(paste(names(x), x, sep="", collapse=""))
    })
    
    # molDF <- data.frame(do.call(rbind, molList))
    # mf <- icRanalysis::generate_mf(.data=molDF, c_cname = "C", o_cname = "O", h_cname = "H",
    #                                n_cname = "N", s_cname = "S", p_cname = "P")
    return(mf)
}

compound_res$MF <- parseChemicalFormulaField(compound_res$`CHEMICAL-FORMULA`)
#cpd <- "CPD-18750"

saveRDS(compound_res, file.path(data.dir, "compounds.RDS"))
