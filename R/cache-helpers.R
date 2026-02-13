#' @importFrom dplyr %>%
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom zen4R ZenodoManager ZenodoRecord
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames

utils::globalVariables(c(
    "doi"
))


## creates new ZenodoManager object 

.onLoad <- function(libname, pkgname) {
    zenodom <<- ZenodoManager$new(
        token = "bUrBXBEHjZvWHjQ0LKejAhc3Ex8d1utXDUZf4JrOYrhMpykNTcNuP0RzUihd"
    )

}

## creates file cache on user's machine
makeCache <- function() {
    cache_file_path <- '~/AnnotatedBCGEData'
    if(!file.exists(cache_file_path)) {
        bfc <- BiocFileCache(cache_file_path, ask=FALSE)
        return()
    } else {
        return()
    }
}


## downloads file from Zenodo via zen4R package
downloadZenFile <- function(conceptDOI, filepath, zen, v=NULL) {
    if (!is.null(v)) {
        record <- zen$getRecordByConceptDOI(conceptDOI)
        versions <- record$getVersions()
        recDOIs <- versions %>%
            pull(doi)
        
        if (v>length(recDOIs)) {
            stop(paste0("Either the requested version of the data does not",
                        " exist, or there was a syntax error. Try passing the version",
                        " as an int."))
        }
        
        newDOI <- recDOIs[v]
        new_rec <- zen$getRecordByDOI(newDOI)
        new_rec$downloadFiles(path=filepath, quiet=TRUE)
        return()
        
    }
    
    record <- zen$getRecordByConceptDOI(conceptDOI)
    record$downloadFiles(path=filepath, quiet=TRUE)
    return()
}


## makes SE object
seConstructor <- function(exp_path, meta_path) {
    exp_data <- read_tsv(exp_path) %>%
        filterRepeatRows()
    start_col <- colnames(exp_data)[7]
    end_col <- colnames(exp_data)[ncol(exp_data)]
    expression_matrix = makeDataMatrix(exp_data, start_col, end_col)
    
    metadata <- read_tsv(meta_path) %>%
        column_to_rownames('Sample_ID')
    
    feature_data <- makeFeatureData(exp_data)
    
    se = makeSummarizedExperiment(expression_matrix, feature_data, metadata)
    return(se)
}


## helper called for downloading specified version of Zenodo file
chooseVersion <- function(datasetID, v, zen, version, conceptDOI) {
    dir_path <- paste0('~/AnnotatedBCGEData/', datasetID, 'v', v)
    if (!dir.exists(dir_path)) {
        if (v<version) {
            dir.create(dir_path)
            message(paste0("Either the data was not found in the cache or a ",
                           "different version was requested. Downloading now."))
            downloadZenFile(conceptDOI, dir_path, zen, v)
            
            
        } else {
            stop(paste0("By default, the most recent version of the data is ",
                        "downloaded. You do not need to specify a version for ",
                        "the latest version."))
        }
    }
    
    exp_filepath <- paste0('~/AnnotatedBCGEData/', datasetID, 'v', v, '/', datasetID, '.tsv.gz')
    meta_filepath <- paste0('~/AnnotatedBCGEData/', datasetID, 'v', v, '/', datasetID, '_metadata.tsv')
    se <- seConstructor(exp_filepath, meta_filepath)
    return(se)
}


## function for checking most recent version of zen file
checkVersions <- function(conceptDOI, zen) {
    record <- zen$getRecordByConceptDOI(conceptDOI)
    versions <- record$getVersions() %>%
        pull(version)
    max_vers <- max(versions)
    return(max_vers)
    
}