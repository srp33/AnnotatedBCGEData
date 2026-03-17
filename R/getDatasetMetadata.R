#' @importFrom dplyr %>% rename select filter

utils::globalVariables(c(
    "accession_id",
    "publishing_platform",
    "experiment_type",
    "summary",
    "title",
    "overall_design")
)


#' A function that takes a vector of dataset IDs and returns a tibble with 
#'  information about how the data was collected and the overall study.
#'  
#' @param datasets is a vector of strings of the dataset IDs to search for.
#' @return Tibble with information about the experiment that produced the data.
#' @examples getDatasetMetadata(c("GSE41197", "GSE10797"))
#' @export 
getDatasetMetadata <- function(datasets) {
    dataset_meta_path <- file.path(tempdir(), 'dataset_meta.tsv')
    if (!file.exists(dataset_meta_path)) {
        downloadZenFile('10.5281/zenodo.17780657', tempdir())
    }
    dataset_meta <- read_tsv(dataset_meta_path) %>%
        select(accession_id, 
               publishing_platform, 
               experiment_type, 
               title, 
               summary, 
               overall_design) %>%
        rename(Dataset_ID=accession_id) %>%
        rename(Source=publishing_platform) %>%
        rename(Experiment_Type=experiment_type) %>%
        rename(Title=title) %>%
        rename(Summary=summary) %>%
        rename(Overall_Design=overall_design) %>%
        filter(Dataset_ID %in% datasets)
    
    return(dataset_meta)
}