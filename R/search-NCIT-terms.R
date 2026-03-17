#' @importFrom stringr regex str_detect
#' @importFrom dplyr rename filter select pull %>% bind_rows
#' @importFrom tidyr separate_rows

utils::globalVariables(c(
    "code",
    "Preferred Label",
    "Synonyms",
    "URI",
    "Definitions",
    "label"
))

NCIT_identifiers <- c("10.5281/zenodo.17488901", 
                      "NCIT_definitions_filtered.tsv.gz")


#' Function that searches NCIT definitions file. 
#' 
#' Takes a term to search by and the type of the term. This could be the name of
#'  the term itself, the ontology's definition, the associated code, or the URI.
#'  The function searches across all terms and returns matches
#'   to the search with the term name, URI, code, and definition. 
#' 
#' @param String to search by.
#' @param Ontology field to search: Name, URI, Code, or Definition.
#' @return A tibble with a row for each matching ontology term.
#' @examples searchOntologyTerms("Progesterone Receptor Status")
#' @export 
searchOntologyTerms <- function(term, term_type="Name") {
    acceptable_terms <- c("Name", "Definition", "URI", "Code")
    
    if (!(term_type %in% acceptable_terms)) {
        stop(paste0(term_type, 
                    " is not an acceptable search term.",
                    " Try: ", paste(acceptable_terms, collapse=", ")))
    }
    
    def_path <- file.path(tempdir(), NCIT_identifiers[2])
    if (!file.exists(def_path)) {
        downloadZenFile(NCIT_identifiers[1], tempdir())
    }
    NCIT_defs <- read_tsv(def_path)
    
    if (term_type == "Name") {
        df <- searchNames(term, NCIT_defs)
        return(df)
    } else if (term_type == "URI") {
        df <- searchURIs(term, NCIT_defs) 
        return(df)
    } else if (term_type == "Code") {
        df <- searchCodes(term, NCIT_defs)
        return(df)
    } else {
        df <- searchDefinition(term, NCIT_defs)
        return(df)
    }
}


## helper function 1: searchNames, searches through Names in definitions df

searchNames <- function(term, df) {
    ncit_names <- rename(df, label = `Preferred Label`) %>%
        select(label, code)
    
    ncit_synonyms <- rename(df, label = Synonyms) %>%
        select(label, code) %>%
        separate_rows(label, sep = "\\|")
    
    names_syns <- bind_rows(ncit_names, ncit_synonyms) %>%
        filter(str_detect(label, regex(term, ignore_case=TRUE))) %>%
        pull(code)
    
    codes <- unique(names_syns)
    
    df_searched <- filter(df, code %in% codes)
    
    return(df_searched)
}


## helper function 2: searchURI, searches through URIs in definitions df

searchURIs <- function(term, df) {
    df_searched <- filter(df, URI==term) 
    return(df_searched)
}


## helper function 3: searchCodes, searches through codes in definitions df

searchCodes <- function(term, df) {
    df_searched <- filter(df, code==term)
    return(df_searched)
}


## helper function 4: searchDefinitions, searches through defs in definitions df

searchDefinition <- function(term, df) {
    df_searched <- filter(df, 
                          str_detect(Definitions, 
                                     regex(term, ignore_case = TRUE)))
    return(df_searched)
}
