#' @section Overview:
#' AnnotatedBCGEData is an R Bioconductor package that provides access 
#' to more than 100 gene-expression datasets generated using microarray and RNA-
#' sequencing technologies. We collected these datsets from public repositories, 
#' including GEO, ArrayExpress, and TCGA. Where feasible, we obtained the and 
#' reprocessed the raw data using computational pipelines and annotations that 
#' were standardized and more updated than those used to process the data 
#' originally. Additionally, we annotated the sample metadata using oncology
#' terms from the National Cancer Institute Thesaurus. This was necessary 
#' because different data submitters used different terminology to describe 
#' breast-cancer patients and samples. Mapping the metadata to ontology terms 
#' makes it possible for users to search the data using these standardized terms. 
#' 
#' - [getDataset()] allows users to get a SummarizedExperiment object of the chosen dataset.
#' - [searchOntologyTerms()] allows users to browse our custom ontology.
#' - [searchForDatasetsByField()] is one function that finds datasets related to a given term.
#' 
#' Read the vignette for more information and detailed examples of the package functions.
#' 
#' @section Examples:
#' 
#' ```
#' # Searching the ontology by a specific term name
#' searchOntologyTerms("Progesterone Receptor Status")
#' 
#' # Find associated datasets in the ontology mappings by a field code
#' searchForDatasetsByField("C16149")
#' 
#' # Retrieving a dataset
#' getDataset("GSE41197")
#' ```
#' 
#' @keywords internal
"_PACKAGE"


## usethis namespace: start
#' @import rlang
#' @import cli
## usethis namespace: end
NULL 