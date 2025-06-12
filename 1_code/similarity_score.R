library(mapa)

load("2_data/0_object/up_enriched_pathways.rda")

# function for calculating similarity score

pathway_similarity_score <- function(database1, database2, api_key){
  database1_df <- 
    database1 %>% 
    dplyr::mutate(
      pathway_id        = pathway_id,
      pathway_name      = pathway_name,
      describtion       = description,
      pathway_class     = pathway_class,
      p_value           = 0,
      p_adjust          = 0,
      all_id            = str_replace_all(compound_list, fixed("{}"), ";"),
      all_number        = ifelse(is.na(all_id), 0, str_count(all_id, ";") + 1),
      mapped_id         = all_id,
      mapped_number     = all_number,
      mapped_percentage = 1
    ) %>% 
    dplyr::select(
      pathway_id, 
      pathway_name, 
      describtion, 
      pathway_class, 
      p_value, 
      p_adjust, 
      all_id, 
      all_number, 
      mapped_id, 
      mapped_number, 
      mapped_percentage
    )
  
  database2_df <- 
    database2 %>% 
    dplyr::mutate(
      pathway_id        = pathway_id,
      pathway_name      = pathway_name,
      describtion       = description,
      pathway_class     = pathway_class,
      p_value           = 0,
      p_adjust          = 0,
      all_id            = str_replace_all(compound_list, fixed("{}"), ";"),
      all_number        = ifelse(is.na(all_id), 0, str_count(all_id, ";") + 1),
      mapped_id         = all_id,
      mapped_number     = all_number,
      mapped_percentage = 1
    ) %>% 
    dplyr::select(
      pathway_id, 
      pathway_name, 
      describtion, 
      pathway_class, 
      p_value, 
      p_adjust, 
      all_id, 
      all_number, 
      mapped_id, 
      mapped_number, 
      mapped_percentage
    )
  
  enriched_pathways@enrichment_hmdb_result@result <- database1_df
  enriched_pathways@enrichment_metkegg_result@result <- database2_df
  
  embeding_res <- 
    mapa::get_bioembedsim(
      object = enriched_pathways,
      api_provider = "openai",
      text_embedding_model = "text-embedding-3-small",
      api_key = api_key, # Insert your api key
      database = c("hmdb", "metkegg"),
      count.cutoff.hmdb = -1,    # Minimum metabolite count cutoff for HMDB pathways
      count.cutoff.metkegg = -1 # Minimum metabolite count cutoff for KEGG pathways
    )
  
  return(embeding_res$sim_matrix)
  
}


