library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)

source("1_code/similarity_score.R")

# Load KEGG and HMDB pathways
kegg_hsa_pathway <- metpath::kegg_hsa_pathway
hmdb_pathway <- metpath::hmdb_pathway

# Convert KEGG pathway to data frame
n <- length(kegg_hsa_pathway@pathway_id)

# Convert description and pathway class to character vectors
describtion_vec <- vapply(
  kegg_hsa_pathway@describtion,
  function(x) {
    ## element might be null
    if (length(x) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      as.character(x[[1]]) # get the first element as character
    }
  },
  FUN.VALUE = character(1)
)

# Extract pathway class, handling potential null elements
pathway_class_vec <- vapply(
  kegg_hsa_pathway@pathway_class,
  function(x) {
    ## element might be null
    if (length(x) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      as.character(x[[1]]) # get the first element as character
    }
  },
  FUN.VALUE = character(1)
)

# Extract pathway class, handling potential null elements
pathway_class_vec <- vapply(
  kegg_hsa_pathway@pathway_class,
  function(x) {
    ## element might be null
    if (length(x) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      as.character(x[[1]]) # get the first element as character
    }
  },
  FUN.VALUE = character(1)
)

compound_list_vec <- vapply(
  kegg_hsa_pathway@compound_list,
  function(x) {
    if (length(x$KEGG.ID) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      paste(x$KEGG.ID, collapse = "{}") # join compound IDs with semicolon
    }
  },
  FUN.VALUE = character(1)
)


# Extract and organize the pathway information into a data frame
kegg_database <- data.frame(
  pathway_id    = kegg_hsa_pathway@pathway_id,
  pathway_name  = kegg_hsa_pathway@pathway_name,
  description   = describtion_vec,
  compound_list = compound_list_vec,
  pathway_class = pathway_class_vec,
  database       = "KEGG"
)

# Replace semicolon with "{}" in pathway_class
kegg_database <- 
  kegg_database %>% 
  dplyr::mutate(
    pathway_class = str_replace_all(pathway_class, "; ", "{}")
  )

save(kegg_database, file = "2_data/7_kegg_database/kegg_database.rda")

# Convert HMDB pathway to data frame
n <- length(hmdb_pathway@pathway_id)

# Convert description and pathway class to character vectors
describtion_vec <- vapply(
  hmdb_pathway@describtion,
  function(x) {
    ## element might be null
    if (length(x) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      as.character(x[[1]]) # get the first element as character
    }
  },
  FUN.VALUE = character(1)
)

pathway_class_vec <- vapply(
  hmdb_pathway@pathway_class,
  function(x) {
    ## element might be null
    if (length(x) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      as.character(x[[1]]) # get the first element as character
    }
  },
  FUN.VALUE = character(1)
)

compound_list_vec <- vapply(
  hmdb_pathway@compound_list,
  function(x) {
    if (length(x$HMDB.ID) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      paste(x$HMDB.ID, collapse = "{}") # join compound IDs with semicolon
    }
  },
  FUN.VALUE = character(1)
)


# Extract and organize the pathway information into a data frame
hmdb_database <- data.frame(
  pathway_id    = hmdb_pathway@pathway_id,
  pathway_name  = hmdb_pathway@pathway_name,
  description   = describtion_vec,
  compound_list = compound_list_vec,
  pathway_class = pathway_class_vec,
  database       = "HMDB"
)

# Replace semicolon with "{}" in pathway_class
hmdb_database <- 
  hmdb_database %>% 
  dplyr::mutate(
    pathway_class = str_replace_all(pathway_class, ";", "{}")
  )

save(hmdb_database, file = "2_data/8_hmdb_database/hmdb_database.rda")

# for test
hmdb_database <- hmdb_database[1:100, ]

similarity_score <- pathway_similarity_score(
  database1 = kegg_database,
  database2 = hmdb_database,
  api_key = "" # Insert your api key
)


# convert the similarity score matrix
smp_rows  <- grepl("^SMP", rownames(similarity_score))   
kegg_cols <- grepl("^hsa", colnames(similarity_score))   

cross_mat <- similarity_score[smp_rows, kegg_cols, drop = FALSE]  

cross_df <- cross_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("SMP_pathway") %>%         
  pivot_longer(-SMP_pathway,                     
               names_to  = "KEGG_pathway",
               values_to = "similarity_score") %>% 
  arrange(SMP_pathway, KEGG_pathway)             

kegg_smpdb_similarity <- cross_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("SMP_pathway") %>%              
  pivot_longer(-SMP_pathway,                         
               names_to  = "KEGG_pathway",
               values_to = "similarity_score") %>% 
  select(KEGG_pathway, SMP_pathway, similarity_score) %>%  
  arrange(KEGG_pathway, SMP_pathway)   

save(kegg_smpdb_similarity, file = "2_data/9_kegg_smpdb_similarity/kegg_smpdb_similarity.rda")








