#' @importFrom zen4R ZenodoManager
NULL

#'  Function that creates a SummarizedExperiment object
#'  
#'  Takes a name of a data set and an optional version of the data, passes
#'  the identifier vector to get concept DOI for specified data set.
#'  
#' @param datasetID the name of the data set, passed as a string
#' @param identifier list of identifiers, defaults to list included
#' @param v version of data to download, default to most recent
#' @return a SummarizedExperiment object for the data set
#' @examples bcgeData("GSE41197")
#' @export
bcgeData <- function(datasetID, identifier=identifiers, zen=zenodom, v=NULL) {
    makeCache()
    
    identifier_vec <- identifier[[datasetID]]
    
    version <- checkVersions(identifier_vec[1], zen)
    
    if (!is.null(v)) {
        conceptDOI <- identifier_vec[1]
        se <- (chooseVersion(datasetID, v, zen, version, conceptDOI))
        return(se)
    }
    
    dir_filepath <- paste0('~/AnnotatedBCGEData/',datasetID, 'v', version)
    exp_filepath <- paste0('~/AnnotatedBCGEData/', datasetID, 'v', version, '/', identifier_vec[2])
    
    if (!dir.exists(dir_filepath)) {
        dir.create(paste0('~/AnnotatedBCGEData/',datasetID, 'v', version))
    }
    
    if (!file.exists(exp_filepath)) {
        conceptDOI <- identifier_vec[1]
        message(paste0("Either the data was not found in the cache or a newer",
                       " version of the data was identified. Downloading now."))
        downloadZenFile(conceptDOI, dir_filepath, zen)
    }
    
    meta_filepath <- paste0('~/AnnotatedBCGEData/', datasetID, 'v', version, '/', identifier_vec[3])
    se <- seConstructor(exp_filepath, meta_filepath)
    
    return(se)
}
