#' Function that creates downloads data and creates SummarizedExperiment object for user to access
#' 
#' Takes a vector with the identifiers for the download URL in the form c(expression data, metadata)
#' 
#' @param identifiers vector of identifiers found in identifiers
#' @return a SummarizedExperiment object for the data set
#' @export
makeObject <- function(identifiers) {
    expression_data = downloadZenodoFile(identifiers) %>%
        filterRepeatRows()
    start_col = colnames(expression_data)[7]
    end_col = colnames(expression_data)[ncol(expression_data)]
    
    expression_matrix = makeDataMatrix(expression_data, start_col, end_col)
    
    sample_metadata = getMetadata(identifiers)
    
    feature_data = makeFeatureData(expression_data)
    
    se = makeSummarizedExperiment(expression_matrix, feature_data, sample_metadata)
    
    return(se)
}


#' named list with the names as the data set accessions and values as a vector with c(expressions, metadata)
#' @format a named list with data sets as names and a vector with URL identifiers as value
#' \describe{
#'  \item{GSE41197}{c(data identifiers, metadata identifiers)}
#'  \item{GSE10797}{c(data identifiers, metadata identifiers)}
#'  \item{GSE59772}{c(data identifiers, metadata identifiers)}
#' }
#' @source generated internally for AnnotatedBCGEData
"identifiers"
identifiers <- list(
    GSE41197 = c('17428998','GSE41197.tsv.gz', '17429158', "GSE41197_metadata.tsv"),
    GSE10797 = c('17429390','GSE10797.tsv.gz', '17429390', 'GSE10797_metadata.tsv'),
    GSE59772 = c('17429395','GSE59772.tsv.gz', '17429395', 'GSE59772_metadata.tsv')
)
