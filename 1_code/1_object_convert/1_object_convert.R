library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(metpath)
library(mapa)

# Load KEGG and HMDB pathways
kegg_hsa_pathway <- metpath::kegg_hsa_pathway
hmdb_pathway <- metpath::hmdb_pathway

# Load original object
load("2_data/0_object/up_enriched_pathways.rda")

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

# Extract all KEGG IDs and their counts from the compound list
all_id <- vapply(
  kegg_hsa_pathway@compound_list,
  function(x) {
    if (length(x$KEGG.ID) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      paste(x$KEGG.ID, collapse = ";") # join compound IDs with semicolon
    }
  },
  FUN.VALUE = character(1)
)

all_number <- vapply(
  kegg_hsa_pathway@compound_list,
  function(x) {
    if (length(x$KEGG.ID) == 0) {
      0L                   # if the element is empty, return 0
    } else {
      length(x$KEGG.ID)   # count the number of KEGG IDs
    }
  },
  FUN.VALUE = integer(1)
)

# Extract and organize the pathway information into a data frame
kegg_has_df <- data.frame(
  pathway_id    = kegg_hsa_pathway@pathway_id,
  pathway_name  = kegg_hsa_pathway@pathway_name,
  describtion   = describtion_vec,
  pathway_class = pathway_class_vec,
  p_value       = rep(0, n),
  p_adjust      = rep(0, n),
  all_id        = all_id,
  all_number    = all_number,
  mapped_id     = all_id,
  mapped_number = all_number,
  mapped_percentage = rep(1, n)
  )

save(kegg_has_df, file = "2_data/1_kegg_hsa/kegg_hsa_pathway.rda")

# Convert hmdb_pathway to data frame
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

# Extract all HMDB IDs and their counts from the compound list
all_id <- vapply(
  hmdb_pathway@compound_list,
  function(x) {
    if (length(x$HMDB.ID) == 0) {
      NA_character_        # if the element is empty, return NA
    } else {
      paste(x$HMDB.ID, collapse = ";") # join compound IDs with semicolon
    }
  },
  FUN.VALUE = character(1)
)

all_number <- vapply(
  hmdb_pathway@compound_list,
  function(x) {
    if (length(x$HMDB.ID) == 0) {
      0L                   # if the element is empty, return 0
    } else {
      length(x$HMDB.ID)   # count the number of HMDB IDs
    }
  },
  FUN.VALUE = integer(1)
)

# Extract and organize the pathway information into a data frame
hmdb_df <- data.frame(
  pathway_id    = hmdb_pathway@pathway_id,
  pathway_name  = hmdb_pathway@pathway_name,
  describtion   = describtion_vec,
  pathway_class = pathway_class_vec,
  p_value       = rep(0, n),
  p_adjust      = rep(0, n),
  all_id        = all_id,
  all_number    = all_number,
  mapped_id     = all_id,
  mapped_number = all_number,
  mapped_percentage = rep(1, n)
)

save(hmdb_df, file = "2_data/2_hmdb/hmdb_pathway.rda")

# Revise the original object
enriched_pathways@enrichment_metkegg_result@result <- kegg_has_df
enriched_pathways@enrichment_hmdb_result@result <- hmdb_df

# Save the revised object
save(enriched_pathways, file = "2_data/3_revised_object/object.rda")
