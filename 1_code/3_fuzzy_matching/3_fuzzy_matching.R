library(r4projects)
setwd(get_project_wd())
rm(list = ls())

# laod KEGG and HMDB pathways
load("2_data/1_kegg_hsa/kegg_hsa_pathway.rda")
load("2_data/2_hmdb/hmdb_pathway.rda")

if (!requireNamespace("fuzzyjoin", quietly = TRUE)) {
  install.packages("fuzzyjoin")      # 其中也会自动安装 stringdist
}

library(dplyr)
library(stringr)
library(fuzzyjoin)

# preprocessing to increase the matching rate
kegg_clean <- 
  kegg_has_df %>% 
  mutate(name_clean = str_remove_all(str_to_lower(pathway_name), "[^a-z0-9]"))

hmdb_clean <-
  hmdb_df %>% 
  mutate(name_clean = str_remove_all(str_to_lower(pathway_name), "[^a-z0-9]"))

# fuzzy matching
threshold <- 0.8 # set a threshold for matching（Jaro-Winkler ≥ 0.85）
max_dist  <- 1 - threshold # set a maximum distance for matching

matched_tbl <- stringdist_inner_join(
  kegg_clean, hmdb_clean,
  by         = c("name_clean" = "name_clean"),  
  method     = "jw",                            # Jaro-Winkler
  max_dist   = max_dist,                       
  distance_col = "jw_dist"                    
) %>% 
  select(kegg_pathway_name = pathway_name.x,
         kegg_pathway_id   = pathway_id.x,
         hmdb_pathway_name = pathway_name.y,
         hmdb_pathway_id   = pathway_id.y,
         jw_dist) %>%                          
  arrange(jw_dist)

# save
save(matched_tbl, file = "2_data/6_fuzzy_matching/matched_pathways.rda")







