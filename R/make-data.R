#' @importFrom zen4R ZenodoManager
NULL

#'  Function that creates a SummarizedExperiment object with breast cancer gene
#'  expression data.
#'  
#'  Takes a name of the desired data set as a string, an optional file path for
#'   the data to be saved to, and an optional version number. The function 
#'   looks for the data either in the specified location or in a temporary 
#'   directory. If the data is not found, it is downloaded. Returns the desired
#'   data set in the form of a SummarizedExperiment object. 
#'  
#' @param datasetID the name of the data set, passed as a string.
#' @param cacheDirPath the optional file path to where the cache should be saved.
#' @param identifier list of identifiers, defaults to list included.
#' @param v version of data to download, default to most recent.
#' @return a SummarizedExperiment object for the data set.
#' @export
bcgeData <- function(datasetID, cacheDirPath=tempdir(), identifier=identifiers, v=NULL) {
    if (cacheDirPath != tempdir()) {
        makeCache(cacheDirPath)
    }
    zen = ZenodoManager$new() 
    
    identifier_vec <- identifier[[datasetID]]
    
    version <- checkVersions(identifier_vec[1])
    
    if (!is.null(v)) {
        conceptDOI <- identifier_vec[1]
        se <- (chooseVersion(datasetID, cacheDirPath, v, version, conceptDOI))
        return(se)
    }
    
    dir_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', version)
    exp_filepath <- paste0(cacheDirPath, '/', 
                           datasetID, 'v', version, '/', identifier_vec[2])
    
    if (!dir.exists(dir_filepath)) {
        dir.create(paste0(cacheDirPath, '/', datasetID, 'v', version))
    }
    
    if (!file.exists(exp_filepath)) {
        conceptDOI <- identifier_vec[1]
        message(paste0("Either the data was not found in the cache or a newer",
                       " version of the data was identified. Downloading now."))
        downloadZenFile(conceptDOI, dir_filepath)
    }
    
    meta_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', 
                            version, '/', identifier_vec[3])
    se <- seConstructor(exp_filepath, meta_filepath)
    
    return(se)
}