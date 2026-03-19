#' Function that searches breast cancer ontology fields mapped to data sets. 
#' 
#' Takes a search code and searches across datasets to find those that are 
#' mapped to that term. 
#' Returns a data frame that shows the name of the data set, the term, and its 
#' code. Users can then pass the names of the datasets to getDataset to 
#' get the SummarizedExperiment objects. 
#' 
#' @param code is the NCIT code that was found from searchOntologyTerms.
#' @return a data frame from the searched file with data set metadata.
#' @examples searchForDatasetFields("C16149")
#' @export
searchForDatasetFields <- function(code) {
    filtered_path <- file.path(tempdir(), 'filtered_mapped_data.tsv.gz')
    
    if (!file.exists(filtered_path)) {
        downloadZenFile('10.5281/zenodo.17583904', tempdir())
    }
    
    filtered_data <- read_tsv(filtered_path)
    df <- searchFields(code, filtered_data)
    return(df)
    
}




