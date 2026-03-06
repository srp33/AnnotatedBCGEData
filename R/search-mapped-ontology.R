#' @importFrom dplyr rename %>% left_join filter

utils::globalVariables(c(
    "accession_id",
    "NCIT_field_code",
    "NCIT_value_code"
))

#' Function that searches breast cancer ontology mapped to data sets. 
#' 
#' Takes a search code and whether to search terms (columns of a particular data
#'  set) or values (values in columns of particular data set). Searches a data 
#' frame that contains data sets, the term associated with it, the values 
#' associated with it, and the terms and values as they were originally found in
#'  the actual data. 
#' Returns a data frame joined with the data set metadata. This shows the name 
#' of the data set, the terms and values, and some metadata about where the data
#'  came from. Users can then pass the names of the data sets to getBCGEData to 
#' get the SummarizedExperiment objects. 
#' 
#' @param code NCIT code that was found from searchNCITerms.
#' @param term_type "field" (the column) or "value" (the data).
#' @return a data frame from the searched file with data set metadata.
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