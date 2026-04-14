get_top_level_gene_counts <- function() {
  
  top_level_ids <- unique(reactome_pathway_hierarchy$top_level_pathway)
  top_level_ids <- top_level_ids[!is.na(top_level_ids)]
  
  gene_counts <- data.frame(
    Top_Level_ID   = character(),
    Top_Level_Name = character(),
    Gene_Count     = integer(),
    stringsAsFactors = FALSE
  )
  
  for (id in top_level_ids) {
    
    # Get human-readable name
    top_name <- legend_reactome_feature_set$Pathway_name[
      legend_reactome_feature_set$Pathway_id == id
    ]
    top_name <- ifelse(length(top_name) > 0 && !is.na(top_name[1]), top_name[1], id)
    
    tryCatch({
      participants <- getParticipants(id)  
      
      # Extract unique gene identifiers from the nested refEntities list,
      # keeping only ReferenceGeneProduct (excludes small molecules, ChEBI, etc.)
      all_genes <- unique(unlist(lapply(participants$refEntities, function(ref) {
        ref$identifier[ref$schemaClass %in% c("ReferenceGeneProduct", "ReferenceIsoform")]
      })))
      
      gene_counts <- rbind(gene_counts, data.frame(
        Top_Level_ID   = id,
        Top_Level_Name = top_name,
        Gene_Count     = length(all_genes),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      message(paste("Failed for:", id, "-", e$message))
      gene_counts <<- rbind(gene_counts, data.frame(
        Top_Level_ID   = id,
        Top_Level_Name = top_name,
        Gene_Count     = NA_integer_,
        stringsAsFactors = FALSE
      ))
    })
  }
  
  gene_counts <- gene_counts[order(gene_counts$Gene_Count, decreasing = TRUE), ]
  return(gene_counts)
}

top_level_gene_counts <- get_top_level_gene_counts()
print(top_level_gene_counts)





plot_categories_reactome <- function(list_sig_reactome_pathways_ids, filename = NULL) {
  
  plot_data <- data.frame(Factor = character(), Top_Level_Pathway = character())
  
  for (i in seq_along(list_sig_reactome_pathways_ids)) {
    
    sig_pathways_ids <- list_sig_reactome_pathways_ids[[i]]
    top_names <- character(length(sig_pathways_ids))
    
    for (j in seq_along(sig_pathways_ids)) {
      id <- sig_pathways_ids[j]
      
      if (id %in% reactome_pathway_hierarchy$top_level_pathway) {
        top_id <- id
      } else {
        top_id <- reactome_pathway_hierarchy$top_level_pathway[match(id, reactome_pathway_hierarchy$pathway)]
        
        if (is.na(top_id) || length(top_id) == 0) {
          tryCatch({
            top_id <- getPathways(id, top.level = TRUE)$stId[1]
          }, error = function(e) {
            message(paste("ID not found in file or API:", id))
            top_id <- NA
          })
        }
      }
      
      if (!is.na(top_id)) {
        top_name <- legend_reactome_feature_set$Pathway_name[legend_reactome_feature_set$Pathway_id == top_id]
        top_names[j] <- ifelse(length(top_name) > 0, top_name, NA)
      } else {
        top_names[j] <- NA
      }
    }
    
    plot_data <- rbind(
      plot_data, 
      data.frame(Factor = rep(names(list_sig_reactome_pathways_ids)[i], length(top_names)), 
                 Top_Level_Pathway = top_names)
    )
  }
  
  plot_counts <- plot_data %>%
    filter(!is.na(Top_Level_Pathway)) %>%
    group_by(Factor, Top_Level_Pathway) %>%
    summarise(Count = n(), .groups = "drop")
  
  # --- NEW: append gene count to each top-level pathway name in the legend ---
  plot_counts <- plot_counts %>%
    mutate(Top_Level_Pathway = {
      gene_n <- top_level_gene_counts$Gene_Count[
        match(Top_Level_Pathway, top_level_gene_counts$Top_Level_Name)
      ]
      ifelse(!is.na(gene_n),
             paste0(Top_Level_Pathway, " (", gene_n, " genes)"),
             Top_Level_Pathway)
    })
  # --------------------------------------------------------------------------
  
  p <- ggplot(plot_counts, aes(x = Factor, y = Count, fill = Top_Level_Pathway)) +
    geom_bar(stat = "identity", position = "fill") +  
    labs(title = "", 
         x = "", 
         y = "", 
         fill = "Top-Level Pathway") +
    scale_x_discrete(labels = names(list_sig_reactome_pathways_ids)) +
    my_theme
  
  if (!is.null(filename)) 
    ggsave(filename = paste0("ADD/", filename, ".pdf"), plot = p, width = 20, height = 7, dpi = 300)
  
  return(p)
}