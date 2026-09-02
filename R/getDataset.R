#' @importFrom zen4R ZenodoManager
NULL

#' Function that retrieves a breast cancer gene expression dataset.

#' Accepts the identifier of a desired dataset, an optional directory path for
#'   the data to be cached in, and an optional version number. The function 
#'   looks for the data in the cache. If the data is not found, it is downloaded
#'   and cached. The data are packaged as SummarizedExperiment objects. If the
#'   user specifies a cache directory, downloaded files are stored in this
#'   location so they don't need to be re-downloaded. The default is to use the
#'   R session's temporary directory, which is removed after the session ends.
#'  
#' @param datasetID The identifier of the data set.
#' @param cacheDirPath The optional directory path to where the cache should be
#'  saved.
#' @param v Version of data to download. Defaults to most recent.
#' @return A SummarizedExperiment object for the dataset.
#' @examples getDataset("GSE41197")
#' @export
getDataset <- function(datasetID, cacheDirPath=tempdir(), v=NULL) {
    if (cacheDirPath != tempdir()) {
        makeCache(cacheDirPath)
    }
    
    identifier_vec <- getIdentifiers(datasetID)
    
    version <- checkVersions(identifier_vec[2])
    
    if (!is.null(v)) {
        conceptDOI <- identifier_vec[2]
        se <- (chooseVersion(datasetID, cacheDirPath, v, version, conceptDOI))
        return(se)
    }
    
    dir_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', version)
    exp_filepath <- paste0(cacheDirPath, '/', 
                           datasetID, 'v', version, '/', identifier_vec[3])
    
    if (!dir.exists(dir_filepath)) {
        dir.create(paste0(cacheDirPath, '/', datasetID, 'v', version))
    }
    
    if (!file.exists(exp_filepath)) {
        conceptDOI <- identifier_vec[2]
        message("Either the data was not found in the cache or a newer",
                       " version of the data was identified. Downloading now.")
        downloadZenFile(conceptDOI, dir_filepath)
    }
    
    meta_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', 
                            version, '/', identifier_vec[4])
    se <- seConstructor(exp_filepath, meta_filepath)
    
    return(se)
}