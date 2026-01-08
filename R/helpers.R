#' @importFrom dplyr filter %>% select distinct
#' @importFrom utils download.file
#' @importFrom stringr str_c str_starts
#' @importFrom readr read_tsv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom tibble column_to_rownames
#' @import tidyverse

utils::globalVariables(c(
    'Chromosome',
    'Dataset_ID',
    'Ensembl_Gene_ID',
    'Entrez_Gene_ID',
    'Gene_Biotype',
    'HGNC_Symbol'
))


##  download file from Zenodo
downloadZenodoFile <- function(identifier) {
    tmp_file_path <- str_c(tempdir(), "/", identifier[2])
    
    if (!file.exists(tmp_file_path)) {
        download.file(paste0(
            "https://zenodo.org/records/", 
            identifier[1], 
            "/files/", 
            identifier[2], 
            "?download=1"), 
            tmp_file_path, 
            mode='wb')
    }
    
    return(read_tsv(tmp_file_path))
}


##  filter out repeat rows
filterRepeatRows <- function(expression_matrix) {
    expression_matrix <- expression_matrix %>%
        filter(!str_starts(Chromosome, 'H'))
    
    return(expression_matrix)
}

##  get metadata
getMetadata <- function(identifier) {
    metadata_identifier <- identifier[3:4]
    sample_metadata <- downloadZenodoFile(metadata_identifier) %>%
        column_to_rownames('Sample_ID')
    
    return(sample_metadata)
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
    GSE41197 = c('17428998','GSE41197.tsv.gz', 
                 '17429158', "GSE41197_metadata.tsv"),
    GSE10797 = c('17429390','GSE10797.tsv.gz', 
                 '17429390', 'GSE10797_metadata.tsv'),
    GSE59772 = c('17429395','GSE59772.tsv.gz', 
                 '17429395', 'GSE59772_metadata.tsv'),
    GSE10281 = c('17665426', 'GSE10281.tsv.gz',
                 '17665426', 'GSE10281_metadata.tsv'),
    GSE10780 = c("17665492", "GSE10780.tsv.gz",
                 "17665492", "GSE10780_metadata.tsv"),
    GSE96058_NextSeq = c("17665296", "GSE96058_NextSeq.tsv.gz",
                         "17665296", "GSE96058_NextSeq_metadata.tsv"),
    GSE96058_HiSeq = c("", "GSE96058_HiSeq.tsv.gz",
                       "", "GSE96058_HiSeq_metadata.tsv"),
    GSE10810 = c('17858011', 'GSE10810.tsv.gz',
                 '17858011', 'GSE10810_metadata.tsv'),
    GSE11001 = c('17858150', 'GSE11001.tsv.gz',
                 '17858150','GSE11001_metadata.tsv'),
    GSE11121 = c('17858268', 'GSE11121.tsv.gz',
                 '17858268', 'GSE11121_metadata.tsv'),
    GSE9574 = c('17869995', 'GSE9574.tsv.gz',
                '17869995', 'GSE9574_metadata.tsv'),
    GSE111662 = c('17869962', 'GSE111662.tsv.gz',
                  '17869962', 'GSE111662_metadata.tsv'),
    GSE93332 = c('17870368', 'GSE93332.tsv.gz',
                 '17870368', 'GSE93332_metadata.tsv'),
    GSE118432 = c('17870380', 'GSE118432.tsv.gz',
                  '17870380', "GSE118432_metadata.tsv"),
    GSE120129 = c('17870745', 'GSE120129.tsv.gz',
                  '17870745', 'GSE120129_metadata.tsv'),
    GSE12093 = c('17870941', 'GSE12093.tsv.gz',
                 '17870941', 'GSE12093_metadata.tsv'),
    GSE12276 = c('17871034', 'GSE12276.tsv.gz',
                 '17871034', 'GSE12276_metadata.tsv'),
    GSE12763 = c("17871140", 'GSE12763.tsv.gz',
                 '17871140', 'GSE12763_metadata.tsv'),
    GSE13787 = c('17871261', 'GSE13787.tsv.gz',
                 '17871261', 'GSE13787_metadata.tsv'),
    GSE9195 = c('17870562', 'GSE9195.tsv.gz',
                '17870562', 'GSE9195_metadata.tsv'),
    GSE14017 = c('17871467', 'GSE14017.tsv.gz',
                 '17871467', 'GSE14017_metadata.tsv'),
    GSE90521 = c('17871540', 'GSE90521.tsv.gz',
                 '17871540', 'GSE90521_metadata.tsv'),
    GSE14018 = c('17871561', 'GSE14018.tsv.gz',
                 '17871561', 'GSE14018_metadata.tsv'),
    GSE8977 = c('17871598', 'GSE8977.tsv.gz',
                '17871598', 'GSE8977_metadata.tsv'),
    GSE1456_U133A = c('17871632', 'GSE1456_U133A.tsv.gz',
                      '17871632', 'GSE1456_U133A_metadata.tsv'),
    GSE1456_U133B = c('17871707', 'GSE1456_U133B.tsv.gz',
                      '17871707', 'GSE1456_U133B_metadata.tsv'),
    GSE86374 = c('17871705', 'GSE86374.tsv.gz',
                 '17871705', 'GSE86374_metadata.tsv'),
    GSE1561 = c('17871757', 'GSE1561.tsv.gz',
                '17871757', 'GSE1561_metadata.tsv'),
    GSE8193 = c('17871782', 'GSE8193.tsv.gz',
                '17871782', 'GSE8193_metadata.tsv'),
    GSE16391 = c('17871811', 'GSE16391.tsv.gz',
                 '17871811', 'GSE16391_metadata.tsv'),
    GSE16446 = c('18126026', 'GSE16446.tsv.gz',
                 '18126026', 'GSE16446_metadata.tsv'),
    GSE167213 = c('18125981', 'GSE167213.tsv.gz',
                  '18125981', 'GSE167213_metadata.tsv.gz'),
    GSE16873 = c('18125929', 'GSE16873.tsv.gz',
                 '18125929', 'GSE16873_metadata.tsv'),
    GSE17705 = c('18125873', 'GSE17705.tsv.gz',
                 '18125873', 'GSE17705_metadata.tsv')
)

