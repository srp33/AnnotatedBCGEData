# code that creates metadata.csv in inst/extdata
# columns = title, description, BiocVersion, Genome, SourceType, SourceURL, SourceVersion, Species, TaxonomyID, 
#   Coordinate_1_based, DataProvider, Maintainer

# 107 data sets total

# REQUIRED WHEN A NEW DATA SET IS ADDED:
#   - INCREASE num_datasets
#   - ADD GEO ACCESSION TO TITLE
#   - ADD GEO NAME TO DESCRIPTION
#   - ADD OSF URL TO SOURCE URL

library(tidyverse)

num_datasets = 3
meta <- data.frame(
    Title = c(paste0('GSE41197'),
              paste0('GSE59772'),
              paste0('GSE10797')
    ),
    Description = c(
        paste0('Differentially Expressed Genes Regulating the Progression ',
               'of Ductal Carcinoma In Situ to Invasive Breast Cancer'),
        paste0('Analysis of compartment-specific gene expression in breast cancer tumors'),
        paste0('Transcriptomes of breast epithelium and stroma in normal reduction ', 
               'mammoplasty and invasive breast cancer patients.')
    ),
    BiocVersion = rep('3.2',num_datasets),
    SourceType = rep('tsv.gz',num_datasets),
    SourceURL = c(
        paste0('https://osf.io/pdc8h'),
        paste0('https://osf.io/ps2kb'),
        paste0('https://osf.io/ebycg')
    ),
    SourceVersion = rep('April 22, 2025',num_datasets),
    Species = rep('Homo Sapiens', num_datasets),
    TaxonomyID = rep(9606, num_datasets),
    DataProvider = rep('OSF', num_datasets),
    Maintainer = rep('Bioconductor Package Maintainer <maintainer@bioconductor.org>', num_datasets)
)

write.csv(meta, file='inst/extdata/metadata.csv', row.names = FALSE)