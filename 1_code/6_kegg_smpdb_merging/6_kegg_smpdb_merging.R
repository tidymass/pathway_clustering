library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

# Load the data
load("2_data/7_kegg_database/kegg_database.rda")
load("2_data/8_hmdb_database/hmdb_database.rda")
load("2_data/9_kegg_smpdb_similarity/kegg_smpdb_similarity.rda")

# for test
hmdb_database <- hmdb_database[1:100, ]

# set the threshold cutoff
threshold <- 0.6

# filter and keep the highest similarity_score per KEGG_pathway
matched_similarity <- kegg_smpdb_similarity %>%
  filter(similarity_score >= threshold) %>%
  group_by(KEGG_pathway) %>%
  slice_max(order_by = similarity_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# merge matched kegg and smpdb with id
merged_df <- 
  matched_similarity %>%
  left_join(kegg_database, by = c("KEGG_pathway" = "pathway_id")) %>%
  left_join(hmdb_database, by = c("SMP_pathway" = "pathway_id"), suffix = c("_1", "_2"))

# combine columns
kegg_smpdb_matched <- 
  merged_df %>%
  rowwise() %>%
  mutate(
    pathway_id = paste(KEGG_pathway, SMP_pathway, sep = "{}"),
    pathway_name = paste(
      unique(na.omit(c(
        str_split(pathway_name_1, fixed("{}"))[[1]],
        str_split(pathway_name_2, fixed("{}"))[[1]]
      ))),
      collapse = "{}"
    ),
    description = paste(
      unique(na.omit(c(
        str_split(description_1, fixed("{}"))[[1]],
        str_split(description_2, fixed("{}"))[[1]]
      ))),
      collapse = "{}"
    ),
    compound_list = paste(
      unique(na.omit(c(
        str_split(compound_list_1, fixed("{}"))[[1]],
        str_split(compound_list_2, fixed("{}"))[[1]]
      ))),
      collapse = "{}"
    ),
    pathway_class = paste(
      unique(na.omit(c(
        str_split(pathway_class_1, fixed("{}"))[[1]],
        str_split(pathway_class_2, fixed("{}"))[[1]]
      ))),
      collapse = "{}"
    ),
    database = "KEGG{}SMPDB"
    
  ) %>% 
  ungroup() %>% 
  select(
    pathway_id, pathway_name, description, compound_list, pathway_class, database
  )

save(kegg_smpdb_matched, file = "2_data/10_kegg_smpdb_merge/kegg_smpdb_matched.rda")

#



# obtain new pathway_name and description
kegg_smpdb_matched <- llm_interpret(
  database = kegg_smpdb_matched,
  api_key = "" # Insert your api key
)




