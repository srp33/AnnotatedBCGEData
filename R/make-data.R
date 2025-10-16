library(tidyverse)
library(SummarizedExperiment)

#' Function that creates downloads data and creates SummarizedExperiment object for user to access
#' 
#' Takes a vector with the identifiers for the download URL in the form c(expression data, metadata)
#' 
#' @param identifiers vector of identifiers found in identifier_list
#' @return a SummarizedExperiment object for the data set
#' @export
makeObject = function(identifiers) {
    data_identifier = identifiers[1]
    meta_identifier = identifiers[2]
    
    expression_data = downloadOSFFile(data_identifier) %>%
        filterRepeatRows()
    start_col = colnames(expression_data)[7]
    end_col = colnames(expression_data)[ncol(expression_data)]
    
    expression_matrix = makeDataMatrix(expression_data, start_col, end_col)
    
    sample_metadata = getMetadata(meta_identifier)
    
    feature_data = makeFeatureData(expression_data)
    
    se = makeSummarizedExperiment(expression_matrix, feature_data, sample_metadata)
    
    return(se)
}


#' named list with the names as the data set accessions and values as a vector with c(expressions, metadata)
#' @format a named list with data sets as names and a vector with URL identifiers as value
#' \describe{
#'  \item{GSE41197}{c(data identifier, metadata identifier)}
#'  \item{GSE10797}{c(data identifier, metadata identifier)}
#'  \item{GSE59772}{c(data identifier, metadata identifier)}
#' }
#' @source generated internally for AnnotatedBCGEData
"identifiers"
identifiers <- list(
    GSE41197 = c('pdc8h','u7x9k'),
    GSE10797 = c('ebycg','vmhuj'),
    GSE59772 = c('ps2kb','dc3qh')
)