library(SummarizedExperiment)
library(tidyverse)

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


#download file from OSF
downloadOSFFile = function(identifier, out_file_path) {
    tmp_file_path = str_c(tempdir(), "/", identifier, ".tsv.gz")
    
    if (!file.exists(tmp_file_path)) {
        download.file(paste0("https://osf.io/download/", identifier), tmp_file_path, mode='wb')
    }
    
    return(read_tsv(tmp_file_path))
}


#filter out repeat rows
filterRepeatRows = function(expression_matrix) {
    expression_matrix = expression_matrix %>%
        filter(!str_starts(Chromosome, 'H'))
    
    return(expression_matrix)
}

# get metadata
getMetadata = function(identifier) {
    sample_metadata = downloadOSFFile(identifier) %>%
        column_to_rownames('Sample_ID')
    
    return(sample_metadata)
}

# make feature data
makeFeatureData = function(expression_matrix) {
    feature_data = select(expression_matrix, Dataset_ID, Entrez_Gene_ID, 
                          HGNC_Symbol, Ensembl_Gene_ID, Chromosome, Gene_Biotype) %>%
        distinct(Ensembl_Gene_ID, .keep_all=TRUE) %>%
        column_to_rownames('Ensembl_Gene_ID')
    
    return(feature_data)
}


# creating expression data matrix
makeDataMatrix = function(dataset, start_col, end_col) {
    expressions = select(dataset, Ensembl_Gene_ID, start_col:end_col)
    expression_matrix = expressions %>%
        column_to_rownames('Ensembl_Gene_ID') %>%
        as.matrix()
    
    return(expression_matrix)
}

#build SummarizedExperiment
makeSummarizedExperiment = function(expressions, features, meta) {
    se = SummarizedExperiment(
        assays = list(counts=expressions),
        rowData = features,
        colData = meta
    )
    
    return(se)
}