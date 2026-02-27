#' @importFrom dplyr rename %>% left_join filter

utils::globalVariables(c(
    "accession_id",
    "NCIT_field_code",
    "NCIT_value_code"
))

#' Function that searches mapped ontology data
#' 
#' Takes a search code and column to search by
#' Returns data frame joined with dataset metadata
#' 
#' @param code NCIT code that was found from searching defs
#' @param term_type field or value
#' @param zen arg for ZenodoManager object
#' @return a data frame from the searched file with dataset metadata
#' @export
searchForDatasets <- function(code, term_type="field") {
    zen <- ZenodoManager$new()
    acceptable_terms <- c("field", "value", "Field", "Value", "FIELD", "VALUE")
    
    if (!(term_type %in% acceptable_terms)) {
        stop(paste0(term_type, 
                    " is not an acceptable term.",
                    " Try 'field' or 'value'"))
    }
    
    filtered_path <- file.path(tempdir(), 'filtered_mapped_data.tsv.gz')
    dataset_meta_path <- file.path(tempdir(), 'dataset_meta.tsv')
    
    if (!file.exists(filtered_path)) {
        downloadZenFile('10.5281/zenodo.17583904', tempdir())
    }
    if (!file.exists(dataset_meta_path)) {
        downloadZenFile('10.5281/zenodo.17780657', tempdir())
    }
    
    filtered_data <- read_tsv(filtered_path)
    dataset_meta <- read_tsv(dataset_meta_path)
    
    if (term_type %in% c("field", "Field", "FIELD")) {
        df <- searchFields(code, filtered_data)
        
    } else {
        df <- searchValues(code, filtered_data)
        
    }
    
    dataset_meta <- rename(dataset_meta, dataset = accession_id)
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