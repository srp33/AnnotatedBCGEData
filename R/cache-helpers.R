#' @importFrom dplyr %>% filter select distinct
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom zen4R ZenodoManager ZenodoRecord
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames
#' @importFrom rappdirs user_config_dir
#' @importFrom rstudioapi selectDirectory
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
    'HGNC_Symbol'
))


## creates new ZenodoManager object 
.onLoad <- function(libname, pkgname) {
    suppressMessages({
        zenodom <<- zen4R::ZenodoManager$new(
            token = 
                "bUrBXBEHjZvWHjQ0LKejAhc3Ex8d1utXDUZf4JrOYrhMpykNTcNuP0RzUihd"
        )
    })
    

}


##asks the user where to save the cache
cacheDir <- function() {
    message(paste0("Data is stored in a cache to prevent repeated downloads.",
                   "A location needs to be selected for the cache directory."))
    if (requireNamespace("rstudioapi", quietly=TRUE) &&
        rstudioapi::isAvailable()) {
        dir <- rstudioapi::selectDirectory()
    } else {
        dir <- readline(prompt="Enter directory path: ")
        check_valid <- dir.create(dir, recursive=TRUE, showWarnings=FALSE)
        while (!check_valid) {
            message(paste0("The directory ", dir, " is invalid. Try again."))
            dir <- readline(prompt="Enter a directory path: ")
            check_valid <- dir.create(dir, recursive=TRUE, showWarnings=FALSE)
        }
    }
    
    config_dir <- file.path(rappdirs::user_config_dir("AnnotatedBCGEData"),
                            "config.txt")
    dir.create(dirname(config_dir), recursive=TRUE, showWarnings=FALSE)
    writeLines(dir, config_dir)
    
    return(dir)
}

loadCache <- function() {
    config_file <- file.path(rappdirs::user_config_dir("AnnotatedBCGEData"), 
                             "config.txt")
    if (file.exists(config_file)) {
        loc <- readLines(config_file)
        cache_path <- paste0(loc, '/AnnotatedBCGEData')
        return(cache_path)
    } else {
        return(NULL)
    }
}


## creates file cache on user's machine
makeCache <- function() {
    if (interactive()) {
        config_file <- file.path(rappdirs::user_config_dir("AnnotatedBCGEData"),
                                 "config.txt")
        if (file.exists(config_file)) {
            cache_file_path <- loadCache()
        } else {
            cache_file_path <- paste0(cacheDir(), '/AnnotatedBCGEData')
        }
    } else {
        cache_file_path <- tempdir()
    }
    
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
                        " exist, or there was a syntax error. Try passing the",
                        " version as an int."))
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
    expression_matrix <- makeDataMatrix(exp_data, start_col, end_col)
    
    metadata <- read_tsv(meta_path) %>%
        column_to_rownames('Sample_ID')
    
    feature_data <- makeFeatureData(exp_data)
    
    se <- makeSummarizedExperiment(expression_matrix, feature_data, metadata)
    return(se)
}


## helper called for downloading specified version of Zenodo file
chooseVersion <- function(datasetID, v, zen, version, conceptDOI) {
    if (interactive()) {
        cache_dir <- loadCache()    
    } else {
        cache_dir <- tempdir()
    }
    
    dir_path <- paste0(cache_dir, '/', datasetID, 'v', v)
    if (!dir.exists(dir_path)) {
        dir.create(dir_path)
        message(paste0("Either the data was not found in the cache or a ",
                    "different version was requested. Downloading now."))
        downloadZenFile(conceptDOI, dir_path, zen, v)
    }
    
    exp_filepath <- paste0(cache_dir, '/', datasetID, 'v', v, '/', 
                           datasetID, '.tsv.gz')
    meta_filepath <- paste0(cache_dir, '/', datasetID, 'v', v, '/', 
                            datasetID, '_metadata.tsv')
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

## make list of identifiers
identifiers <- list(
    GSE41197 = c('10.5281/zenodo.17428997','GSE41197.tsv.gz', 
                 "GSE41197_metadata.tsv"),
    GSE10797 = c('10.5281/zenodo.17429389','GSE10797.tsv.gz', 
                 'GSE10797_metadata.tsv'),
    GSE59772 = c('10.5281/zenodo.17429394','GSE59772.tsv.gz', 
                 'GSE59772_metadata.tsv'),
    GSE10281 = c('10.5281/zenodo.17665283', 'GSE10281.tsv.gz',
                 'GSE10281_metadata.tsv'),
    GSE10780 = c("10.5281/zenodo.17665491", "GSE10780.tsv.gz",
                 "GSE10780_metadata.tsv"),
    GSE96058_NextSeq = c("10.5281/zenodo.17665295", "GSE96058_NextSeq.tsv.gz",
                         "GSE96058_NextSeq_metadata.tsv"),
    GSE96058_HiSeq = c("10.5281/zenodo.17665546", "GSE96058_HiSeq.tsv.gz",
                       "GSE96058_HiSeq_metadata.tsv"),
    GSE10810 = c('10.5281/zenodo.17858010', 'GSE10810.tsv.gz',
                 'GSE10810_metadata.tsv'),
    GSE11001 = c('10.5281/zenodo.17858149', 'GSE11001.tsv.gz',
                 'GSE11001_metadata.tsv'),
    GSE11121 = c('10.5281/zenodo.17858267', 'GSE11121.tsv.gz',
                 'GSE11121_metadata.tsv'),
    GSE9574 = c('10.5281/zenodo.17869994', 'GSE9574.tsv.gz',
                'GSE9574_metadata.tsv'),
    GSE111662 = c('10.5281/zenodo.17869961', 'GSE111662.tsv.gz',
                  'GSE111662_metadata.tsv'),
    GSE93332 = c('10.5281/zenodo.17870367', 'GSE93332.tsv.gz',
                 'GSE93332_metadata.tsv'),
    GSE118432 = c('10.5281/zenodo.17870379', 'GSE118432.tsv.gz',
                  "GSE118432_metadata.tsv"),
    GSE120129 = c('10.5281/zenodo.17870744', 'GSE120129.tsv.gz',
                  'GSE120129_metadata.tsv'),
    GSE12093 = c('10.5281/zenodo.17870940', 'GSE12093.tsv.gz',
                 'GSE12093_metadata.tsv'),
    GSE12276 = c('10.5281/zenodo.17871033', 'GSE12276.tsv.gz',
                 'GSE12276_metadata.tsv'),
    GSE12763 = c("10.5281/zenodo.17871139", 'GSE12763.tsv.gz',
                 'GSE12763_metadata.tsv'),
    GSE13787 = c('10.5281/zenodo.17871260', 'GSE13787.tsv.gz',
                 'GSE13787_metadata.tsv'),
    GSE9195 = c('10.5281/zenodo.17870561', 'GSE9195.tsv.gz',
                'GSE9195_metadata.tsv'),
    GSE14017 = c('10.5281/zenodo.17871466', 'GSE14017.tsv.gz',
                 'GSE14017_metadata.tsv'),
    GSE90521 = c('10.5281/zenodo.17871539', 'GSE90521.tsv.gz',
                 'GSE90521_metadata.tsv'),
    GSE14018 = c('10.5281/zenodo.17871560', 'GSE14018.tsv.gz',
                 'GSE14018_metadata.tsv'),
    GSE8977 = c('10.5281/zenodo.17871597', 'GSE8977.tsv.gz',
                'GSE8977_metadata.tsv'),
    GSE1456_U133A = c('10.5281/zenodo.17871631', 'GSE1456_U133A.tsv.gz',
                      'GSE1456_U133A_metadata.tsv'),
    GSE1456_U133B = c('10.5281/zenodo.17871706', 'GSE1456_U133B.tsv.gz',
                      'GSE1456_U133B_metadata.tsv'),
    GSE86374 = c('10.5281/zenodo.17871704', 'GSE86374.tsv.gz',
                 'GSE86374_metadata.tsv'),
    GSE1561 = c('10.5281/zenodo.17871756', 'GSE1561.tsv.gz',
                'GSE1561_metadata.tsv'),
    GSE8193 = c('10.5281/zenodo.17871781', 'GSE8193.tsv.gz',
                'GSE8193_metadata.tsv'),
    GSE16391 = c('10.5281/zenodo.17871810', 'GSE16391.tsv.gz',
                 'GSE16391_metadata.tsv'),
    GSE16446 = c('10.5281/zenodo.18126025', 'GSE16446.tsv.gz',
                 'GSE16446_metadata.tsv'),
    GSE167213 = c('10.5281/zenodo.18125980', 'GSE167213.tsv.gz',
                  'GSE167213_metadata.tsv.gz'),
    GSE16873 = c('10.5281/zenodo.18125928', 'GSE16873.tsv.gz',
                 'GSE16873_metadata.tsv'),
    GSE17705 = c('10.5281/zenodo.18125872', 'GSE17705.tsv.gz',
                 'GSE17705_metadata.tsv'),
    GSE17907 = c('10.5281/zenodo.18125795', 'GSE17907.tsv.gz',
                 'GSE17907_metadata.tsv'),
    GSE18728 = c('10.5281/zenodo.18125752', 'GSE18728.tsv.gz',
                 'GSE18728_metadata.tsv'),
    GSE18864 = c('10.5281/zenodo.18125678', 'GSE18864.tsv.gz',
                 'GSE18864_metadata.tsv'),
    GSE19615 = c('10.5281/zenodo.18125616', 'GSE19615.tsv.gz',
                 'GSE19615_metadata.tsv'),
    GSE19697 = c('10.5281/zenodo.18125500', 'GSE19697.tsv.gz',
                 'GSE19697_metadata.tsv'),
    GSE20086 = c('10.5281/zenodo.18125411', 'GSE20086.tsv.gz',
                 'GSE20086_metadata.tsv'),
    GSE20181 = c('10.5281/zenodo.18125340', 'GSE20181.tsv.gz',
                 'GSE20181_metadata.tsv'),
    GSE20194 = c('10.5281/zenodo.18125265', 'GSE20194.tsv.gz',
                 'GSE20194_metadata.tsv'),
    GSE20271 = c('10.5281/zenodo.18125215', 'GSE20271.tsv.gz',
                 'GSE20271_metadata.tsv'),
    GSE2034 = c('10.5281/zenodo.18125180', 'GSE2034.tsv.gz',
                'GSE2034_metadata.tsv'),
    GSE20437 = c('10.5281/zenodo.18125148', 'GSE20437.tsv.gz',
                 'GSE20437_metadata.tsv'),
    GSE20685 = c('10.5281/zenodo.18125056', 'GSE20685.tsv.gz',
                 'GSE20685_metadata.tsv'),
    GSE20711 = c('10.5281/zenodo.18124993', 'GSE20711.tsv.gz',
                 'GSE20711_metadata.tsv'),
    GSE21422 = c('10.5281/zenodo.18124934', 'GSE21422.tsv.gz',
                 'GSE21422_metadata.tsv'),
    GSE21653 = c('10.5281/zenodo.18124875', 'GSE21653.tsv.gz',
                 'GSE21653_metadata.tsv'),
    GSE21947 = c('10.5281/zenodo.18124811', 'GSE21947.tsv.gz',
                 'GSE21947_metadata.tsv'),
    GSE22093 = c('10.5281/zenodo.18124065', 'GSE22093.tsv.gz',
                 'GSE22093_metadata.tsv'),
    GSE22513 = c('10.5281/zenodo.18124025', 'GSE22513.tsv.gz',
                 'GSE22513_metadata.tsv'),
    GSE22544 = c('10.5281/zenodo.18123958', 'GSE22544.tsv.gz',
                 'GSE22544_metadata.tsv'),
    GSE23720 = c('10.5281/zenodo.18123873', 'GSE23720.tsv.gz',
                 'GSE23720_metadata.tsv'),
    GSE23988 = c('10.5281/zenodo.18123776', 'GSE23988.tsv.gz',
                 'GSE23988_metadata.tsv'),
    GSE24185 = c('10.5281/zenodo.18123613', 'GSE24185.tsv.gz',
                 'GSE24185_metadata.tsv'),
    GSE25055 = c('10.5281/zenodo.18123419', 'GSE25055.tsv.gz',
                 'GSE25055_metadata.tsv'),
    GSE25065 = c('10.5281/zenodo.18123307', 'GSE25065.tsv.gz',
                 'GSE25065_metadata.tsv'),
    GSE25407 = c('10.5281/zenodo.18123231', 'GSE25407.tsv.gz',
                 'GSE25407_metadata.tsv'),
    GSE2603 = c('10.5281/zenodo.18123167', 'GSE2603.tsv.gz',
                'GSE2603_metadata.tsv'),
    GSE26910 = c('10.5281/zenodo.18122869', 'GSE26910.tsv.gz',
                 'GSE26910_metadata.tsv'),
    GSE27562 = c('10.5281/zenodo.18122800', 'GSE27562.tsv.gz',
                 'GSE27562_metadata.tsv'),
    GSE28821 = c('10.5281/zenodo.18122749', 'GSE28821.tsv.gz',
                 'GSE28821_metadata.tsv'),
    GSE29431 = c('10.5281/zenodo.18122693', 'GSE29431.tsv.gz',
                 'GSE29431_metadata.tsv'),
    GSE2990 = c('10.5281/zenodo.18122626', 'GSE2990.tsv.gz',
                'GSE2990_metadata.tsv'),
    GSE31138 = c('10.5281/zenodo.18122572', 'GSE31138.tsv.gz',
                 'GSE31138_metadata.tsv'),
    GSE31192 = c('10.5281/zenodo.18122512', 'GSE31192.tsv.gz',
                 'GSE31192_metadata.tsv'),
    GSE31448 = c('10.5281/zenodo.18122451', 'GSE31448.tsv.gz',
                 'GSE31448_metadata.tsv'),
    GSE31519 = c('10.5281/zenodo.18122416', 'GSE31519.tsv.gz', 
                 'GSE31519_metadata.tsv'),
    GSE32518 = c('10.5281/zenodo.17971892', 'GSE32518.tsv.gz',
                 'GSE32518_metadata.tsv'),
    GSE32646 = c('10.5281/zenodo.17971814', 'GSE32646.tsv.gz',
                 'GSE32646_metadata.tsv'),
    GSE33692 = c('10.5281/zenodo.17971725', 'GSE33691.tsv.gz',
                 'GSE33692_metadata.tsv'),
    GSE3494_U133A = c('10.5281/zenodo.17971566', 'GSE3494_U133A.tsv.gz',
                      'GSE3494_U133A_metadata.tsv'),
    GSE3494_U133B = c('10.5281/zenodo.17971509', 'GSE3494_U133B.tsv.gz',
                      'GSE3494_U133B_metadata.tsv'),
    GSE41194 = c('10.5281/zenodo.17971457', 'GSE41194.tsv.gz',
                 'GSE41194_metadata.tsv'),
    GSE41196 = c('10.5281/zenodo.17971410', 'GSE41196.tsv.gz',
                 'GSE41196_metadata.tsv'),
    GSE42568 = c('10.5281/zenodo.17971357', 'GSE42568.tsv.gz',
                 'GSE42568_metadata.tsv'),
    GSE43365 = c('10.5281/zenodo.17971290', 'GSE43365.tsv.gz',
                 'GSE43365_metadata.tsv'),
    GSE45255 = c('10.5281/zenodo.17971228', 'GSE45255.tsv.gz',
                 'GSE45255_metadata.tsv'),
    GSE4611 = c('10.5281/zenodo.17971104', 'GSE4611.tsv.gz',
                'GSE4611_metadata.tsv'),
    GSE46184 = c('10.5281/zenodo.17971043', 'GSE46184.tsv.gz',
                 'GSE64184_metadata.tsv'),
    GSE47109 = c('10.5281/zenodo.17970992', 'GSE47109.tsv.gz',
                 'GSE47109_metadata.tsv'),
    GSE48390 = c('10.5281/zenodo.17970913', 'GSE48390.tsv.gz',
                 'GSE48390_metadata.tsv'),
    GSE4922_U133A = c('10.5281/zenodo.17970839', 'GSE4922_U133A.tsv.gz',
                      'GSE4922_U133A_metadata.tsv'),
    GSE4922_U133B = c('10.5281/zenodo.17970768', 'GSE4922_U133B.tsv.gz',
                      'GSE4922_U133B_metadata.tsv'),
    GSE50948 = c('10.5281/zenodo.17917747', 'GSE50948.tsv.gz',
                 'GSE50948_metadata.tsv'),
    GSE5327 = c('10.5281/zenodo.17917717', 'GSE5327.tsv.gz',
                'GSE5327_metadata.tsv'),
    GSE5460 = c('10.5281/zenodo.17917686', 'GSE5460.tsv.gz',
                'GSE5460_metadata.tsv'),
    GSE5764 = c('10.5281/zenodo.17917665', 'GSE5764.tsv.gz',
                'GSE5764_metadata.tsv'),
    GSE57968 = c('10.5281/zenodo.17917638', 'GSE57968.tsv.gz',
                 'GSE57968_metadata.tsv'),
    GSE5847 = c('10.5281/zenodo.17917623', 'GSE5847.tsv.gz',
                'GSE5847_metadata.tsv'),
    GSE58644 = c('10.5281/zenodo.17917594', 'GSE58644.tsv.gz',
                 'GSE58644_metadata.tsv'),
    GSE58984 = c('10.5281/zenodo.17917572', 'GSE58984.tsv.gz',
                 'GSE58984_metadata.tsv'),
    GSE61304 = c('10.5281/zenodo.17917546', 'GSE61304.tsv.gz',
                 'GSE61304_metadata.tsv'),
    GSE62944_Normal = c('10.5281/zenodo.17917528', 'GSE62944_Normal.tsv.gz',
                        'GSE62944_Normal_metadata.tsv'),
    GSE62944_Tumor = c('10.5281/zenodo.17917495', 'GSE62944_Tumor.tsv.gz',
                       'GSE62944_Tumor_metadata.tsv'),
    GSE6434 = c('10.5281/zenodo.17872512', 'GSE6434.tsv.gz',
                'GSE6434_metadata.tsv'),
    GSE6532_U133A = c('10.5281/zenodo.17872456', 'GSE6532_U133A.tsv.gz',
                      'GSE6532_U133A_metadata.tsv'),
    GSE6532_U133B = c('10.5281/zenodo.17872391', 'GSE6532_U133B.tsv.gz',
                      'GSE6532_U133B_metadata.tsv'),
    GSE6532_U133Plus2 = c('10.5281/zenodo.17872340', 'GSE6532_U133Plus2.tsv.gz',
                          'GSE6532_U133Plus2_metadata.tsv'),
    GSE7378 = c('10.5281/zenodo.17872270', 'GSE7378.tsv.gz',
                'GSE7378_metadata.tsv'),
    GSE7390 = c('10.5281/zenodo.17872246', 'GSE7390.tsv.gz',
                'GSE7390_metadata.tsv'),
    GSE76275 = c('10.5281/zenodo.17872182', 'GSE76275.tsv.gz',
                 'GSE76275_metadata.tsv'),
    GSE7904 = c('10.5281/zenodo.17872097', 'GSE7904.tsv.gz',
                'GSE7904_metadata.tsv'),
    GSE81538 = c('10.5281/zenodo.17872045', 'GSE81538.tsv.gz',
                 'GSE81538_metadata.tsv'),
    GSE81838 = c('10.5281/zenodo.17871832', 'GSE81838.tsv.gz',
                 'GSE81838_metadata.tsv'),
    E_TABM_158 = c('10.5281/zenodo.18126090', 'E_TABM_158.tsv.gz',
                   'E_TABM_158_metadata.tsv'),
    ICGC_KR = c('10.5281/zenodo.18272171', 'ICGC_KR.tsv.gz',
                'ICGC_KR.metadata.tsv'),
    METABRIC = c('10.5281/zenodo.18272234', 'METABRIC.tsv.gz',
                 'METABRIC_metadata.tsv')
    
)