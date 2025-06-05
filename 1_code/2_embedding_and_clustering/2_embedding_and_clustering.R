library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(mapa)

# load enriched pathways object
load("2_data/3_revised_object/object.rda")
object <- enriched_pathways

# for test purpose, here select the top 100 hmdb pathways
object@enrichment_hmdb_result@result <- 
  object@enrichment_hmdb_result@result[1:100,]

# merge pathways via embedding
embeding_res <- 
  mapa::get_bioembedsim(
  object = object,
  api_provider = "openai",
  text_embedding_model = "text-embedding-3-small",
  api_key = "", # Insert your api key
  database = c("hmdb", "metkegg"),
  count.cutoff.hmdb = 0,    # Minimum metabolite count cutoff for HMDB pathways
  count.cutoff.metkegg = -1 # Minimum metabolite count cutoff for KEGG pathways
)

save(embeding_res, file = "2_data/4_embedding_results/embeding_res.rda")

# clustering
cluster_res <- mapa::merge_pathways_bioembedsim(
  object = embeding_res,
  cluster_method = "girvan newman",
  sim.cutoff = 0.6 # similarity cutoff for pathways
)

# for test use
object_res <- cluster_res
object_res@merged_module$functional_module_result <- object@merged_module$functional_module_result[1:3,]

# llm interpretation
llm_interpret_res <- mapa::llm_interpret_module(
  object = object,
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = "",
  embedding_output_dir = "/Users/yijiang/Tidymass/Metabolic_pathway/2_data/6_embedding_output"
)


