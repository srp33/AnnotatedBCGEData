#' @importFrom dplyr rename %>% left_join filter

utils::globalVariables(c(
    "geo_accession"
))

filtered_data_url <- "https://zenodo.org/records/17603341/files/filtered_mapped_data.tsv.gz?download=1"
dataset_meta_url <- "https://zenodo.org/records/17780658/files/combined_data.tsv?download=1"

#' Function that searches mapped ontology data
#' 
#' Takes a search code and column to search by
#' Returns data frame joined with dataset metadata
#' 
#' @param code NCIT code that was found from searching defs
#' @param term_type field or value
#' @return a data frame from the searched file with dataset metadata
#' @examples searchForDatasets("C19790", "field")
#' @export
searchForDatasets <- function(code, term_type="field") {
    acceptable_terms <- c("field", "value", "Field", "Value", "FIELD", "VALUE")
    
    if (!(term_type %in% acceptable_terms)) {
        stop(paste0(term_type, 
                    " is not an acceptable term.",
                    " Try 'field' or 'value'"))
    }
    
    filtered_data <- downloadZenodoFile(c("17603341", "filtered_mapped_data.tsv.gz"))
    dataset_meta <- downloadZenodoFile(c("17780658", "combined_data.tsv"))
    
    if (term_type %in% c("field", "Field", "FIELD")) {
        df <- searchFields(code, filtered_data)
        
    } else {
        df <- searchValues(code, filtered_data)
        
    }
    
    dataset_meta <- rename(dataset_meta, dataset = geo_accession)
    df_searched_meta <- left_join(df, dataset_meta)
    
    return(df_searched_meta)
    
}


## helper function that searches just the field codes
searchFields <- function(code, df) {
    df_searched <- filter(df, NCIT_field_code == code)
    return (df_searched)
}

## helper function that searches just the value codes
searchValues <- function(code, df) {
    df_searched <- filter(df, NCIT_value_code == code)
    return (df_searched)
}