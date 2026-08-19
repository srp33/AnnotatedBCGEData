#' @importFrom dplyr %>% filter select distinct pull rename
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom zen4R ZenodoManager ZenodoRecord
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames
#' @importFrom stringr str_starts
#' @importFrom readr read_tsv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom tibble column_to_rownames

utils::globalVariables(c(
    "doi",
    'Chromosome',
    'Dataset_ID',
    'Ensembl_Gene_ID',
    'Entrez_Gene_ID',
    'Gene_Biotype',
    'HGNC_Symbol',
    'accession_id',
    "NCIT_field_code",
    "NCIT_value_code",
    "NCIT_field",
    "NCIT_values",
    "dataset", 
    "orig_field",
    "orig_values",
    "Expressions_DOI",
    "Expressions_Filename",
    "Metadata_Filename"
))


## checks if the filepath exists and a cache can be made
tryPath <- function(path) {
    tryCatch(
        {
            bfc <- BiocFileCache(path, ask=FALSE)
            return(TRUE)
        },
        
        error = function(cond) {
            return(FALSE)
        }
    )
}


## checks if there needs to be sleep time between requests
tryRequestConcept <- function(conceptDOI) {
    tryCatch(
        {
            zen <- ZenodoManager$new()
            rec <- zen$getRecordByConceptDOI(conceptDOI)
            if (!inherits(rec, "ZenodoRecord")) {
                stop()
            }
            return(list(rec, TRUE))
        },
        
        error = function(cond) {
            noRec <- "not a record"
            return(list(noRec, FALSE))
        }
    )
}

tryRequestRecord <- function(recDOI) {
    tryCatch(
        {
            zen <- ZenodoManager$new()
            rec <- zen$getRecordByDOI(recDOI)
            if(!inherits(rec, "ZenodoRecord")) {
                stop()
            }
            return(list(rec, TRUE))
        },
        
        error = function(cond) {
            noRec <- "not a record"
            return(list(noRec, FALSE))
        }
    )
}


# trycatch block for getVersions
tryGetVersions <- function(rec) {
    tryCatch(
        {
            version_attempt <- rec$getVersions()
            if(inherits(version_attempt, "ZenodoException")) {
                stop()
            }
            return(list(version_attempt, TRUE))
        },
        
        error = function(cond) {
            noRec <- "not a version list"
            return(list(noRec, FALSE))
        }
    )
}


## creates file cache on user's machine
makeCache <- function(cacheDirPath) {
    if(!file.exists(cacheDirPath)) {
        if(tryPath(cacheDirPath)) {
            return()  
        } else {
            stop("The filepath is invalid! Are you sure the path exists?")
        }
    } else {
        return()
    }
}


## downloads file from Zenodo via zen4R package
downloadZenFile <- function(conceptDOI, filepath, v=NULL) {
    recCheck <- tryRequestConcept(conceptDOI)
    record <- recCheck[[1]]
    CHECK <- recCheck[[2]]
    sleep_time <- 0
    
    while (CHECK != TRUE) {
        if (sleep_time >= 60) {
            stop("Maximum sleep exceeded")
        }
        Sys.sleep(15)
        message("####################################\n",
                       "Please wait...\n", "Due to Zenodo limitations, we are",
                       " only allowed to submit a certain number of requests ",
                       "per minute. To avoid exceeding those limits, we are ",
                       "pausing for 15 seconds.\n",
                       "####################################")
        sleep_time <- sleep_time + 15
        recCheck <- tryRequestConcept(conceptDOI)
        record <- recCheck[[1]]
        CHECK <- recCheck[[2]]
    }
    
    if (!is.null(v)) {
        version_attempt <- tryGetVersions(record)
        version_CHECK <- version_attempt[[2]]
        versions <- version_attempt[[1]]
        sleep_time <- 0
        
        while (version_CHECK != TRUE) {
            if (sleep_time >= 60) {
                stop("Maximum sleep time exceeded")
            }
            Sys.sleep(15)
            message("####################################\n",
                           "Please wait...\n", "Due to Zenodo limitations, we are",
                           " only allowed to submit a certain number of requests ",
                           "per minute. To avoid exceeding those limits, we are ",
                           "pausing for 15 seconds.\n",
                           "####################################")
            sleep_time <- sleep_time + 15
            version_attempt <- tryGetVersions(record)
            version_CHECK <- version_attempt[[2]]
            versions <- version_attempt[[1]]
        }
        recDOIs <- versions %>%
            dplyr::pull(doi)
        
        if (v>length(recDOIs)) {
            stop("Either the requested version of the data does not",
                        " exist, or there was a syntax error. Try passing the",
                        " version as an int.")
        }
        
        newDOI <- recDOIs[v]
        recAttempt <- tryRequestRecord(newDOI)
        new_rec <- recAttempt[[1]]
        CHECK <- recAttempt[[2]]
        sleep_time <- 0
        
        while (CHECK != TRUE) {
            if (sleep_time >= 60) {
                stop("Maximum sleep exceeded")
            }
            Sys.sleep(15)
            message("####################################\n",
                           "Please wait...\n", "Due to Zenodo limitations, we are",
                           " only allowed to submit a certain number of requests ",
                           "per minute. To avoid exceeding those limits, we are ",
                           "pausing for 15 seconds.\n",
                           "####################################")
            sleep_time <- sleep_time + 15
            recAttempt <- tryRequestRecord(newDOI)
            new_rec <- recAttempt[[1]]
            CHECK <- recAttempt[[2]]
        }
        
        new_rec$downloadFiles(path=filepath, quiet=TRUE)
        return()
        
    }
    
    record$downloadFiles(path=filepath, quiet=TRUE)
    return()
}


## makes SE object
seConstructor <- function(exp_path, meta_path) {
    exp_data <- read_tsv(exp_path) %>%
        filterRepeatRows()
    start_col <- colnames(exp_data)[7]
    end_col <- colnames(exp_data)[ncol(exp_data)]
    expression_matrix <- makeDataMatrix(exp_data, start_col, end_col)
    
    metadata <- read_tsv(meta_path) %>%
        column_to_rownames('Sample_ID')
    
    feature_data <- makeFeatureData(exp_data)
    
    se <- makeSummarizedExperiment(expression_matrix, feature_data, metadata)
    return(se)
}


## helper called for downloading specified version of Zenodo file
chooseVersion <- function(datasetID, cacheDirPath, v, version, conceptDOI) {
    dir_path <- paste0(cacheDirPath, '/', datasetID, 'v', v)
    if (!dir.exists(dir_path)) {
        dir.create(dir_path)
        message("Either the data was not found in the cache or a ",
                    "different version was requested. Downloading now.")
        downloadZenFile(conceptDOI, dir_path, v)
    }
    
    exp_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', v, '/', 
                           datasetID, '.tsv.gz')
    meta_filepath <- paste0(cacheDirPath, '/', datasetID, 'v', v, '/', 
                            datasetID, '_metadata.tsv')
    se <- seConstructor(exp_filepath, meta_filepath)
    return(se)
}


## function for checking most recent version of zen file
checkVersions <- function(conceptDOI) {
    recCheck <- tryRequestConcept(conceptDOI)
    rec <- recCheck[[1]]
    CHECK <- recCheck[[2]]
    sleep_time <- 0
    
    while (CHECK != TRUE) {
        if (sleep_time >= 60) {
            stop("Maximum sleep exceeded")
        }
        Sys.sleep(15)
        message("####################################\n",
                       "Please wait...\n", "Due to Zenodo limitations, we are",
                       " only allowed to submit a certain number of requests ",
                       "per minute. To avoid exceeding those limits, we are ",
                       "pausing for 15 seconds.\n",
                       "####################################")
        sleep_time <- sleep_time + 15
        recCheck <- tryRequestConcept(conceptDOI)
        rec <- recCheck[[1]]
        CHECK <- recCheck[[2]]
    }
    
    version_attempt <- tryGetVersions(rec)
    version_CHECK <- version_attempt[[2]]
    versions <- version_attempt[[1]]
    sleep_time <- 0
    
    while (version_CHECK != TRUE) {
        if (sleep_time >= 60) {
            stop("Maximum sleep time exceeded")
        }
        Sys.sleep(15)
        message("####################################\n",
                       "Please wait...\n", "Due to Zenodo limitations, we are",
                       " only allowed to submit a certain number of requests ",
                       "per minute. To avoid exceeding those limits, we are ",
                       "pausing for 15 seconds.\n",
                       "####################################")
        sleep_time <- sleep_time + 15
        version_attempt <- tryGetVersions(rec)
        version_CHECK <- version_attempt[[2]]
        versions <- version_attempt[[1]]
    }
    
    versions <- versions %>%
        dplyr::pull(version)
    max_vers <- max(versions)
    return(max_vers)
    
}


##  filter out repeat rows
filterRepeatRows <- function(expression_matrix) {
    expression_matrix <- expression_matrix %>%
        filter(!str_starts(Chromosome, 'H'))
    
    return(expression_matrix)
}


##  make feature data
makeFeatureData <- function(expression_matrix) {
    feature_data <- select(expression_matrix, 
                           Dataset_ID, Entrez_Gene_ID, 
                           HGNC_Symbol, Ensembl_Gene_ID, 
                           Chromosome, Gene_Biotype) %>%
        distinct(Ensembl_Gene_ID, .keep_all=TRUE) %>%
        column_to_rownames('Ensembl_Gene_ID')
    
    return(feature_data)
}


##  creating expression data matrix
makeDataMatrix <- function(dataset, start_col, end_col) {
    expressions <- select(dataset, Ensembl_Gene_ID, start_col:end_col)
    expression_matrix <- expressions %>%
        column_to_rownames('Ensembl_Gene_ID') %>%
        as.matrix()
    
    return(expression_matrix)
}

##  build SummarizedExperiment
makeSummarizedExperiment <- function(expressions, features, meta) {
    se <- SummarizedExperiment(
        assays = list(counts=expressions),
        rowData = features,
        colData = meta
    )
    
    return(se)
}

## helper function that searches just the field codes
searchFields <- function(code, df) {
    df_searched <- filter(df, NCIT_field_code %in% code) %>%
        select(dataset, orig_field, NCIT_field_code) %>%
        rename(Dataset_ID=dataset) %>%
        rename(Field=orig_field) %>%
        rename(Code=NCIT_field_code) %>%
        distinct()
    return (df_searched)
}

## helper function that searches just the value codes
searchValues <- function(code, df) {
    df_searched <- filter(df, NCIT_value_code %in% code) %>%
        select(dataset, orig_field, NCIT_field_code, orig_values, NCIT_value_code) %>%
        rename(Dataset_ID=dataset) %>%
        rename(Field=orig_field) %>%
        rename(Field_Code=NCIT_field_code) %>%
        rename(Original_Value=orig_values) %>%
        rename(Code=NCIT_value_code) %>%
        distinct()
    return (df_searched)
}


## find dataset identifiers from list of available datasets
getIdentifiers <- function(datasetID) {
    ident_filepath <- file.path(tempdir(), "identifier_spreadsheet.tsv.txt")
    if (!file.exists(ident_filepath)) {
        downloadZenFile('10.5281/zenodo.21985420', tempdir())
    }
    ident_tib <- read_tsv(ident_filepath) 
    
    if (datasetID %in% ident_tib$Dataset_ID) {
        target_row <- filter(ident_tib, Dataset_ID==datasetID) %>%
            select(Dataset_ID, Expressions_DOI, Expressions_Filename, Metadata_Filename)
        target_vec <- unlist(target_row[1,], use.names=FALSE)
        return(target_vec)
    } else {
        stop("Dataset not found.")
    }
}