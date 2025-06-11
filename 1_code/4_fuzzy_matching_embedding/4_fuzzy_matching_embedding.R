library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(mapa)
library(ggplot2)

load("2_data/6_fuzzy_matching/matched_pathways.rda")
load("2_data/3_revised_object/object.rda")

object <- enriched_pathways

# obtain the most similar pathways 
matched <- 
  matched_tbl[1:40,]

result_new <- 
  result %>% 
  dplyr::filter(pathway_id %in% matched$hmdb_pathway_id)

object@enrichment_hmdb_result@result <- 
  result_new

# merge pathways via embedding
embeding_res <- 
  mapa::get_bioembedsim(
    object = object,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = "", # Insert your api key
    database = c("hmdb", "metkegg"),
    count.cutoff.hmdb = -1,    # Minimum metabolite count cutoff for HMDB pathways
    count.cutoff.metkegg = -1 # Minimum metabolite count cutoff for KEGG pathways
  )

similarity_score <- 
  embeding_res$sim_matrix


similarity_long <- 
  similarity_score %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "hmdb_pathway_id") %>%
  pivot_longer(
    cols = -hmdb_pathway_id,
    names_to = "kegg_pathway_id",
    values_to = "similarity_score"
  )

matched_score <- 
  left_join(matched, similarity_long, by = c("kegg_pathway_id", "hmdb_pathway_id"))

ggplot(matched_score, aes(x = "All Pairs", y = similarity_score)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkblue", size = 1.5) +  # 添加散点
  labs(title = "Similarity Score for Overlapped Pathways in KEGG and SMPDB",
       x = "", y = "Similarity Score") +
  theme_minimal()
