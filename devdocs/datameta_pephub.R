library(tidyverse)
library(yaml)

tar_file_url = "https://pephub.databio.org/geo?view=archive"

tar_file <- "../devdocs/geo_2025_04_28.tar"

extract <- "../devdocs"

#untar(tar_file, exdir=extract)

gse_ids <- c("GSE41197", "GSE10797", "GSE59772", "GSE10281", "GSE10780", 
             "GSE96058", "GSE10810", "GSE11001", "GSE11121", "GSE9574",
             "GSE111662", "GSE93332", "GSE118432", "GSE120129", "GSE12093",
             "GSE12276", "GSE12763", "GSE13787", "GSE9195", "GSE14017",
             "GSE90521", "GSE14018", "GSE8977", "GSE1456", "GSE1561",
             "GSE86374", "GSE8193", "GSE16391", "GSE16446", "GSE167213", 
             "GSE16873", "GSE17705", "GSE17907", "GSE18728", "GSE18864",
             "GSE19615", "GSE19697", "GSE20086", "GSE20181", "GSE20194",
             "GSE20271", "GSE2034", "GSE20437", "GSE20685", "GSE20711",
             "GSE21422", "GSE21653", "GSE21947", "GSE22093", "GSE22513",
             "GSE22544", "GSE23720", "GSE23988", "GSE24185", "GSE25055",
             "GSE25065", "GSE25407", "GSE2603", "GSE26910", "GSE27562",
             "GSE28821", "GSE29431", "GSE2990", "GSE31138", "GSE31192",
             "GSE31448", "GSE31519", "GSE32518", "GSE32646", "GSE33692",
             "GSE3494", "GSE41194", "GSE41196", "GSE42568", "GSE43365",
             "GSE45255", "GSE4611", "GSE46184", "GSE47109", "GSE48390",
             "GSE4922", "GSE50948", "GSE5327", "GSE5460", "GSE5764", 'GSE57968',
             "GSE5847", "GSE58644", "GSE58984", "GSE61304", "GSE62944",
             "GSE6434", "GSE6532", "GSE7378", "GSE7390", "GSE76275", "GSE7904",
             "GSE81538", "GSE81838")

find_values <- function(len, strt, id_vec) {
    matching_ids <- c()
    for (id in id_vec) {
        if (str_length(id)==len) {
            id_nums <- substr(id, 4, len)
            if (str_starts(id_nums, strt)) {
                matching_ids <- c(matching_ids, id)
            }
        }
    }
    return (matching_ids)
}

gse1000s_dir <- "../devdocs/peps/gse1nnn"
target <- "../devdocs/pepsZipUsed"

copy_files <- function(id_vec, target_dir) {
    for (id in id_vec) {
        if (str_length(id)==7){
            id_nums <- substr(id, 4, str_length(id))
            thousand <- substr(id_nums, 1, 1)
            lower_id <- str_to_lower(id)
            source_dir <- paste0("../devdocs/peps/gse", thousand, "nnn/", lower_id, ".zip")
            file.copy(source_dir, target_dir, overwrite=TRUE)
        } else if (str_length(id)==8){
            id_nums <- substr(id, 4, str_length(id))
            thousand <- substr(id_nums, 1,2)
            lower_id <- str_to_lower(id)
            source_dir <- paste0("../devdocs/peps/gse", thousand, "nnn/", lower_id, ".zip")
            file.copy(source_dir, target_dir, overwrite=TRUE)
        } else if (str_length(id)==9) {
            id_nums <- substr(id, 4, str_length(id))
            thousand <- substr(id_nums, 1,3)
            lower_id <- str_to_lower(id)
            source_dir <- paste0("../devdocs/peps/gse", thousand, "nnn/", lower_id, ".zip")
            file.copy(source_dir, target_dir, overwrite=TRUE)
        }
        zipped <- list.files(target_dir)
        for (z in zipped) {
            zip_name <- paste0(target_dir, "/", z)
            unzip(zipfile=zip_name, exdir=target_dir)
        }
        
    }
}

#copy_files(gse_ids, target) 

#zips_to_remove <- list.files(target, pattern="\\.zip$", full.names=TRUE)
#file.remove(zips_to_remove)

meta_tib <- tibble(
    accession_id = character(),
    contact_city = character(),
    contact_country = character(),
    contact_institute = character(),
    available_date = character(),
    last_update_date = character(),
    overall_design = character(),
    summary = character(),
    platform_id = character(),
    publishing_platform = character(),
    title = character(),
    experiment_type = character()
)

read_yaml_files <- function(target_dir, tib) {
    files <- list.files(target_dir)
    for (file in files) {
        config_name <- paste0(target_dir, "/", file)
        config <- yaml.load_file(config_name)
        tib <- add_row(tib,
        accession_id=config$experiment_metadata$series_geo_accession,
        contact_city=config$experiment_metadata$series_contact_city,
        contact_country=config$experiment_metadata$series_contact_country,
        contact_institute=config$experiment_metadata$series_contact_institute,
        available_date=config$experiment_metadata$series_status,
        last_update_date=config$experiment_metadata$series_last_update_date,
        overall_design=config$experiment_metadata$series_overall_design,
        summary=config$experiment_metadata$series_summary,
        platform_id=config$experiment_metadata$series_platform_id,
        publishing_platform="GEO",
        title=config$experiment_metadata$series_title,
        experiment_type=config$experiment_metadata$series_type)
    }
    return(tib)
}

meta_tib <- read_yaml_files(target, meta_tib)

meta_tib <- add_row(meta_tib, 
                    accession_id="GSE19615",
                    contact_city="Manhattan",
                    contact_country="USA",
                    contact_institute="Memorial Sloane Kettering Cancer Center",
                    available_date="Public on January 05, 2010",
                    last_update_date="Jun 06 2022",
                    overall_design=paste0("Affymetrix U133Plus2 from 115 tumor",
                                          " patients"),
                    summary=paste0("Integrated DNA and expression array", 
                    " analysis in primary human breast tumors identified", 
                    " chromosome 8q22 copy number gain and a suite of", 
                    " over-expressed genes in this region highly relevant to ",
                    "subsequent recurrence. 8q copy gain is associated with", 
                    " increased risk of distant metastasis despite adjuvant ",
                    "chemotherapy and multiple genes from 8q22 are ",
                    "overexpressed and have pleiotropic effects, contributing ",
                    "to tumor growth, survival, and reduced killing of tumor",
                    " cell by chemotherapy."),
                    platform_id="GPL570",
                    publishing_platform="GEO",
                    title=paste0("Integrated genomic and function", 
                                 " characterization of 8q22 gain"),
                    experiment_type="Expression profiling by array")
meta_tib <- add_row(meta_tib,
                    accession_id="E-TABM-158",
                    contact_city="Berkeley",
                    contact_country="USA",
                    contact_institute="Lawrence Berkeley National Laboratory",
                    available_date="Public on June 27, 2008",
                    last_update_date="Mar 09 2022",
                    overall_design=paste0("We measured baseline gene expression",
                    "profiles for a set of breast tumors. Note: this experiment",
                    " was modified in June 2008 when the CEL files associated",
                    " with some hybridizations were changed."),
                    summary=paste0("Breast cancers today are of predominantly",
                                   " T1 (0.1 ≥ 2.0 cm) or T2 (>2 ≤ 5 cm)",
                                   " categories due to early diagnosis.",
                                   " Molecular profiling using microarrays has",
                                   " led to the notion of breast cancer as a",
                                   " heterogeneous disease both clinically and",
                                   " molecularly. Given the prognostic power",
                                   " and clinical use of tumor size, the",
                                   " purpose of this study was to search for",
                                   " molecular signatures characterizing",
                                   " clinical T1 and T2. In total 46 samples",
                                   " were included in the discovery dataset.",
                                   " After adjusting for hormone receptor",
                                   " status, lymph node status, grade, and",
                                   " tumor subclass 441 genes were differently",
                                   " expressed between T1 and T2 tumors. Focal",
                                   " adhesion and extracellular matrix ",
                                   "receptor interaction were upregulated in ",
                                   "the smaller tumors while p38MAPK signaling",
                                   " and immune-related pathways were more ",
                                   "dominant in the larger tumors. The T-size",
                                   " signature was then tested on a validation",
                                   " set of 947 breast tumor samples. Using the",
                                   " T-size expression signatures instead of",
                                   " tumor size leads to a significant",
                                   " difference in risk for distant metastases",
                                   " (P < 0.001). If further confirmed, ",
                                   "this molecular signature can be used to ",
                                   "select patients with tumor category T1 who",
                                   " may need more aggressive treatment and ",
                                   "patients with tumor category T2 who may ",
                                   "have less benefit from it."),
                    platform_id="A-AFFY-76",
                    publishing_platform="ArrayExpress",
                    title=paste0("Transcription profiling of human breast ",
                    "cancer samples"))

meta_tib <- add_row(meta_tib,
                    accession_id="METABRIC",
                    contact_city="Cambridge",
                    contact_country="UK",
                    contact_institute=paste0("University of Cambridge,",
                    "Department of Oncology"),
                    available_date="Public on December 21, 2012",
                    last_update_date=NA,
                    overall_design=paste0("DNA and RNA were isolated from ",
                    "samples and hybridized to the Affymetrix SNP 6.0 and ",
                    "Illumina HT-12 v3 platforms for genomic and ",
                    "transcriptional profiling, respectively."),
                    summary=paste0("The elucidation of breast cancer subgroups",
                    " and their molecular drivers requires integrated views of",
                    " the genome and transcriptome from representative numbers",
                    " of patients. We present an integrated analysis of copy ",
                    "number and gene expression in a discovery and validation ",
                    "set of 997 and 995 primary breast tumours, respectively, ",
                    "with long-term clinical follow-up. Inherited variants ",
                    "(copy number variants and single nucleotide polymorphisms)",
                    " and acquired somatic copy number aberrations (CNAs) were",
                    " associated with expression in ~40% of genes, with the ",
                    "landscape dominated by cis- and trans-acting CNAs. By ",
                    "delineating expression outlier genes driven in cis by ",
                    "CNAs, we identified putative cancer genes, including ",
                    "deletions in PPP2R2A, MTAP and MAP2K4. Unsupervised ",
                    "analysis of paired DNA–RNA profiles revealed novel ",
                    "subgroups with distinct clinical outcomes, which ",
                    "reproduced in the validation cohort. These include a ",
                    "high-risk, oestrogen-receptor-positive 11q13/14 cis-acting",
                    "subgroup and a favourable prognosis subgroup devoid of ",
                    "CNAs. Trans-acting aberration hotspots were found to ",
                    "modulate subgroup-specific gene networks, including a ",
                    "TCR deletion-mediated adaptive immune response in the ",
                    "‘CNA-devoid’ subgroup and a basal-specific chromosome 5 ",
                    "deletion-associated mitotic network. Our results provide ",
                    "a novel molecular stratification of the breast cancer ",
                    "population, derived from the impact of somatic CNAs on ",
                    "the transcriptome."),
                    platform_id=NA,
                    publishing_platform="cBioPortal",
                    title="METABRIC")

meta_tib <- add_row(meta_tib,
                    accession_id="ICGC_KR",
                    contact_city=NA,
                    contact_country="Korea",
                    contact_institute= NA,
                    available_date="2012",
                    last_update_date=NA,
                    overall_design=paste0("Breast cancer gene expression data",
                                          " for the International Cancer ",
                                          "Genome Atlas South Korean Cohort"),
                    summary=NA,
                    platform_id=NA,
                    publishing_platform="ICGC",
                    title="ICGC_KR")

write_tsv(meta_tib, "../devdocs/dataset_meta.tsv")