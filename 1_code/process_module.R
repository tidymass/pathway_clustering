# process_module.R

preprocess_module <- function(df,
                              orgdb = org.Hs.eg.db  # Annotation database, can be replaced based on species
) {
  if ("geneID" %in% colnames(df)) {
    # functional_module_result for genes
    query_type <- "gene"
    df <- df %>%
      dplyr::rename(queryid = geneID)
    if ("node" %in% colnames(df)) {
      df <- df %>%
        dplyr::rename(pathway_id = node)
    }
  } else if ("mapped_id" %in% colnames(df)) {
    # functional_module_result for metabolites
    query_type <- "metabolite"
    df <- df %>%
      dplyr::rename(queryid = mapped_id) %>%
      dplyr::rename(pathway_id = node)
  }
  
  required_cols <- c("Description", "queryid", "module", "pathway_id")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Data frame is missing the following columns:", paste(missing_cols, collapse = ", ")))
  }
  
  result_list <- list()
  
  for (i in seq_len(nrow(df))) {
    # Collect module information
    pthName_raw <- df$Description[i]
    qID_raw <- df$queryid[i]
    mod_raw  <- df$module[i]
    pthID_raw <- df$pathway_id[i]
    
    pthName_vec <- unlist(strsplit(pthName_raw, ";"))
    pthID_vec <- unlist(strsplit(pthID_raw, ";"))
    qIDs_vec <- unlist(strsplit(qID_raw, "/"))
    
    pthName_vec <- trimws(pthName_vec)
    qIDs_vec <- trimws(qIDs_vec)
    pthID_vec <- trimws(pthID_vec)
    
    # Retrieve information from databases
    if (query_type == "gene") {
      GeneIDs_vec <- qIDs_vec
      pathway_info <- get_pathway_and_gene_info(pthID_vec)
      pathwayDescription_vec <- pathway_info$pathwayDescription_vec
      pathwayReferencePMID_vec <- pathway_info$PMID_vec
      
      ##这是返回pathway里所有的gene，由于结果太多废弃了
      # GeneIDs_vec <- pathway_info$annotated_entrez_vec
      # GeneSymbols_vec <- pathway_info$annotated_symbol_vec  # 获取标准基因符号
      # GeneDescriptions_vec <- AnnotatedGenename_vec
      ## No need to convert id since id conversion has been completed before this step
      # gene_conversion <- convert_gene_identifiers(query_vec, orgdb)
      # GeneSymbols_vec <- gene_conversion$GeneSymbols_vec
      # GeneNames_vec <- gene_conversion$GeneNames_vec
      
      suppressMessages(GeneSymbols_vec <- AnnotationDbi::mapIds(orgdb, keys = GeneIDs_vec, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first") %>% unname())
      suppressMessages(GeneNames_vec <- AnnotationDbi::mapIds(orgdb, keys = GeneIDs_vec, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first") %>% unname())
      
      result_list[[mod_raw]] <- list(
        PathwayNames = pthName_vec,
        PathwayDescription = pathwayDescription_vec,
        PathwayReferencePMID = unique(pathwayReferencePMID_vec),
        GeneIDs = unique(GeneIDs_vec),
        GeneSymbols = unique(GeneSymbols_vec),
        GeneNames_vec = unique(GeneNames_vec)
      )
      
      # if (mod_raw %in% names(result_list)) {
      #   result_list[[mod_raw]]$PathwayNames <- c(result_list[[mod_raw]]$PathwayNames, pthName_vec)
      #   result_list[[mod_raw]]$PathwayDescription <- c(result_list[[mod_raw]]$PathwayDescription, pathwayDescription_vec)
      #   result_list[[mod_raw]]$PathwayReferencePMID <- unique(c(result_list[[mod_raw]]$PathwayReferencePMID, pathwayReferencePMID_vec))
      #   result_list[[mod_raw]]$GeneIDs <- unique(c(result_list[[mod_raw]]$GeneIDs, GeneIDs_vec))
      #   result_list[[mod_raw]]$GeneSymbols <- unique(c(result_list[[mod_raw]]$GeneSymbols, GeneSymbols_vec))
      #   result_list[[mod_raw]]$GeneNames_vec <- unique(c(result_list[[mod_raw]]$GeneNames_vec, GeneNames_vec))
      # }
    } else if (query_type == "metabolite") {
      MetIDs_vec <- qIDs_vec
      pathway_info <- get_pathway_and_metabolite_info(pthID_vec)
      pathwayDescription_vec = pathway_info$pathwayDescription_vec
      pathwayReferencePMID_vec = pathway_info$PMID_vec
      
      
      MetNames_vec <- get_metabolite_name(MetIDs_vec)
      
      result_list[[mod_raw]] <- list(
        PathwayNames = pthName_vec,
        PathwayDescription = pathwayDescription_vec,
        PathwayReferencePMID = unique(pathwayReferencePMID_vec),
        MetIDs = unique(MetIDs_vec),
        MetNames_vec = unique(MetNames_vec)
      )
    }
  }
  
  return(result_list)
}

# get pathway and metabolite information
get_pathway_and_metabolite_info <- function(pathwayID_vec) {
  # Identify different types of pathway ID
  # kegg_ids <- pathwayID_vec[grepl("^[a-zA-Z]+\\d+$", pathwayID_vec)]
  smpdb_ids <- pathwayID_vec[grepl("^SMP\\d+$", pathwayID_vec)]
  kegg_ids <- pathwayID_vec[!(pathwayID_vec %in% smpdb_ids)]
  
  # Initialize the structure of retrieved data
  pathwayDescription_vec <- c()
  PMID_vec <- c()
  
  # Retrieve KEGG info
  if (length(kegg_ids) > 0) {
    metkegg_info <- get_metkegg_info(kegg_ids)
    for (x in metkegg_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
    }
  }
  
  # Retrieve SMPDB info
  if (length(smpdb_ids) > 0) {
    smpdb_info <- get_smpdb_info(smpdb_ids)
    for (x in smpdb_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
    }
  }
  
  # **返回所有合并后的向量**
  return(list(
    pathwayDescription_vec = pathwayDescription_vec,
    PMID_vec = PMID_vec
  ))
}

# get metabolite name
get_metabolite_name <- function(MetIDs_vec) {
  hmdb_MetNames_vec <-
    metpath::hmdb_compound_database@spectra.info %>%
    dplyr::filter(HMDB.ID %in% MetIDs_vec) %>%
    dplyr::pull(Compound.name)
  
  # kegg_MetNames_vec <-
  #   MetIDs_vec[!grepl("^HMDB", MetIDs_vec)] |>
  #   sapply(function(x){
  #     compound_info <- KEGGREST::keggGet(x)
  #     compound_name <- gsub(";", "", compound_info[[1]]$NAME[1])
  #   })
  kegg_MetNames_vec <-
    MetIDs_vec[!grepl("^HMDB", MetIDs_vec)] |>
    sapply(function(x) {
      # Skip if ID is "NA" (the string)
      if (is.na(x) || x == "NA") {
        return(NA)
      }
      
      # Add delay between requests (1-2 seconds recommended)
      Sys.sleep(1.5)
      
      # Try-catch for error handling
      tryCatch({
        compound_info <- KEGGREST::keggGet(x)
        compound_name <- gsub(";", "", compound_info[[1]]$NAME[1])
        return(compound_name)
      }, error = function(e) {
        cat("Error for ID", x, ":", e$message, "\n")
        return(NA)  # Return NA for failed requests
      })
    })
  
  MetNames_vec <- c(hmdb_MetNames_vec, unname(kegg_MetNames_vec))
  
  # return(MetNames_vec)
  return(unlist(MetNames_vec, use.names = FALSE))
}



get_metkegg_info <- function(kegg_ids) {
  metkegg_info <- list()
  
  for (i in 1:length(kegg_ids)) {
    kegg_id <- kegg_ids[i]
    ### Get a single KEGG entry
    entry <- KEGGREST::keggGet(dbentries = kegg_id)[[1]]
    
    ### Get PMID for reference
    if ("REFERENCE" %in% names(entry)) {
      pmid_list <- lapply(entry$REFERENCE, function(r) {
        if (grepl("^PMID", r$REFERENCE)) {
          pmid <- gsub("[^:]+:", "", r$REFERENCE)
        } else {
          pmid <- NULL
        }
      })
      pmid <- unlist(pmid_list)
    } else {
      pmid <- ""
    }
    
    # Collect base info (always included)
    all_info <- list(
      "id" = unname(entry$ENTRY),
      "term_name" = sub(" - Homo sapiens \\(human\\)$", "", entry$NAME),
      "term_definition" = paste(entry$DESCRIPTION, collapse = " "),
      "PMID" = pmid
    )
    
    metkegg_info <- c(list(all_info), metkegg_info)
  }
  
  return(metkegg_info)
}



get_smpdb_info <- function(smpdb_ids) {
  smpdb_pathway_info <- metpath::hmdb_pathway
  smpdb_info <-
    smpdb_ids %>%
    purrr::map(function(x) {
      indx <- which(smpdb_pathway_info@pathway_id == x)
      term_name <- smpdb_pathway_info@pathway_name[indx]
      term_definition <- smpdb_pathway_info@describtion[indx][[1]]
      annotated_compound_name <- smpdb_pathway_info@compound_list[indx][[1]]$Compound.name
      list(
        "id" = x,
        "term_name" = term_name,
        "term_definition" = term_definition,
        "PMID" = NA_character_
      )
    })
  
  return(smpdb_info)
}
