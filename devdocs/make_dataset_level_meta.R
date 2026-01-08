library(GEOquery)
library(tidyverse)

get_gse_metadata <- function(gse_id) {
    # Get GEO dataset (GSEMatrix = FALSE returns the main metadata structure)
    gse = getGEO(gse_id, GSEMatrix = FALSE)
    
    # Extract metadata list
    meta = Meta(gse)
    
    # Collapse any vector elements into single strings
    meta_collapsed = map(meta, function(x) {
        if (length(x) > 1) paste(x, collapse = "; ") else x
    })
    
    # Create a single-row tibble
    return(tibble::as_tibble_row(meta_collapsed))
}
gse_ids <- c("GSE41197", "GSE10797", "GSE59772", "GSE10281", "GSE10780", 
             "GSE96058", "GSE10810", "GSE11001", "GSE11121", "GSE9574",
             "GSE111662", "GSE93332", "GSE118432", "GSE120129", "GSE12093",
             "GSE12276", "GSE12763", "GSE13787", "GSE9195", "GSE14017",
             "GSE90521", "GSE14018", "GSE8977", "GSE1456", "GSE1561",
             "GSE86374", "GSE8193", "GSE16391", "GSE16446", "GSE167213", 
             "GSE16873", "GSE17705")

all_data <- map(gse_ids, get_gse_metadata)
combined <- bind_rows(all_data) %>%
    select(geo_accession, contact_city, contact_state, contact_country, 
           contact_institute, status, last_update_date, overall_design, 
           summary, platform_id, name, title, type) %>%
    unite(city_state, contact_city, contact_state, sep=", ", na.rm=TRUE) %>%
    unite(study_location, city_state, contact_country, sep=", ") %>%
    rename(available_date = status) %>%
    rename(publishing_platform = name) %>%
    rename(experiment_type = type)



write_tsv(combined, "C:/Users/heidi/Desktop/bioconductor_package/AnnotatedBCGEData/devdocs/combined_data.tsv")
