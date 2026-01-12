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
                 '18125873', 'GSE17705_metadata.tsv'),
    GSE17907 = c('18125796', 'GSE17907.tsv.gz',
                 '18125796', 'GSE17907_metadata.tsv'),
    GSE18728 = c('18125753', 'GSE18728.tsv.gz',
                 '18125753', 'GSE18728_metadata.tsv'),
    GSE18864 = c('18125679', 'GSE18864.tsv.gz',
                 '18125679', 'GSE18864_metadata.tsv'),
    GSE19615 = c('18125617', 'GSE19615.tsv.gz',
                 '18125617', 'GSE19615_metadata.tsv'),
    GSE19697 = c('18125501', 'GSE19697.tsv.gz',
                 '18125501', 'GSE19697_metadata.tsv'),
    GSE20086 = c('18125412', 'GSE20086.tsv.gz',
                 '18125412', 'GSE20086_metadata.tsv'),
    GSE20181 = c('18125341', 'GSE20181.tsv.gz',
                 '18125341', 'GSE20181_metadata.tsv'),
    GSE20194 = c('18125266', 'GSE20194.tsv.gz',
                 '18125266', 'GSE20194_metadata.tsv'),
    GSE20271 = c('18125216', 'GSE20271.tsv.gz',
                 '18125216', 'GSE20271_metadata.tsv'),
    GSE2034 = c('18125181', 'GSE2034.tsv.gz',
                '18125181', 'GSE2034_metadata.tsv'),
    GSE20437 = c('18125149', 'GSE20437.tsv.gz',
                 '18125149', 'GSE20437_metadata.tsv'),
    GSE20685 = c('18125057', 'GSE20685.tsv.gz',
                 '18125057', 'GSE20685_metadata.tsv'),
    GSE20711 = c('18124994', 'GSE20711.tsv.gz',
                 '18124994', 'GSE20711_metadata.tsv'),
    GSE21422 = c('18124935', 'GSE21422.tsv.gz',
                 '18124935', 'GSE21422_metadata.tsv'),
    GSE21653 = c('18124876', 'GSE21653.tsv.gz',
                 '18124876', 'GSE21653_metadata.tsv'),
    GSE21947 = c('18124812', 'GSE21947.tsv.gz',
                 '18124812', 'GSE21947_metadata.tsv'),
    GSE22093 = c('18124066', 'GSE22093.tsv.gz',
                 '18124066', 'GSE22093_metadata.tsv'),
    GSE22513 = c('18124026', 'GSE22513.tsv.gz',
                 '18124026', 'GSE22513_metadata.tsv'),
    GSE22544 = c('18123959', 'GSE22544.tsv.gz',
                 '18123959', 'GSE22544_metadata.tsv'),
    GSE23720 = c('18123874', 'GSE23720.tsv.gz',
                 '18123874', 'GSE23720_metadata.tsv'),
    GSE23988 = c('18123777', 'GSE23988.tsv.gz',
                 '18123777', 'GSE23988_metadata.tsv'),
    GSE24185 = c('18123614', 'GSE24185.tsv.gz',
                 '18123614', 'GSE24185_metadata.tsv'),
    GSE25055 = c('18123420', 'GSE25055.tsv.gz',
                 '18123420', 'GSE25055_metadata.tsv'),
    GSE25065 = c('18123308', 'GSE25065.tsv.gz',
                 '18123308', 'GSE25065_metadata.tsv'),
    GSE25407 = c('18123232', 'GSE25407.tsv.gz',
                 '18123232', 'GSE25407_metadata.tsv'),
    GSE2603 = c('18123168', 'GSE2603.tsv.gz',
                '18123168', 'GSE2603_metadata.tsv'),
    GSE26910 = c('18122870', 'GSE26910.tsv.gz',
                 '18122870', 'GSE26910_metadata.tsv'),
    GSE27562 = c('18122801', 'GSE27562.tsv.gz',
                 '18122801', 'GSE27562_metadata.tsv'),
    GSE28821 = c('18122750', 'GSE28821.tsv.gz',
                 '18122750', 'GSE28821_metadata.tsv'),
    GSE29431 = c('18122694', 'GSE29431.tsv.gz',
                 '18122694', 'GSE29431_metadata.tsv'),
    GSE2990 = c('18122627', 'GSE2990.tsv.gz',
                '18122627', 'GSE2990_metadata.tsv'),
    GSE31138 = c('18122573', 'GSE31138.tsv.gz',
                 '18122573', 'GSE31138_metadata.tsv'),
    GSE31192 = c('18122513', 'GSE31192.tsv.gz',
                 '18122513', 'GSE31192_metadata.tsv'),
    GSE31448 = c('18122452', 'GSE31448.tsv.gz',
                 '18122452', 'GSE31448_metadata.tsv'),
    GSE31519 = c('18122417', 'GSE31519.tsv.gz', 
                 '18122417', 'GSE31519_metadata.tsv'),
    GSE32518 = c('17971893', 'GSE32518.tsv.gz',
                 '17971893', 'GSE32518_metadata.tsv'),
    GSE32646 = c('17971815', 'GSE32646.tsv.gz',
                 '17971815', 'GSE32646_metadata.tsv'),
    GSE33692 = c('17971726', 'GSE33691.tsv.gz',
                 '17971726', 'GSE33692_metadata.tsv'),
    GSE3494_U133A = c('17971567', 'GSE3494_U133A.tsv.gz',
                      '17971567', 'GSE3494_U133A_metadata.tsv'),
    GSE3494_U133B = c('17971510', 'GSE3494_U133B.tsv.gz',
                      '17971510', 'GSE3494_U133B_metadata.tsv'),
    GSE41194 = c('17971458', 'GSE41194.tsv.gz',
                 '17971458', 'GSE41194_metadata.tsv'),
    GSE41196 = c('17971411', 'GSE41196.tsv.gz',
                 '17971411', 'GSE41196_metadata.tsv'),
    GSE42568 = c('17971358', 'GSE42568.tsv.gz',
                 '17971358', 'GSE42568_metadata.tsv'),
    
)

