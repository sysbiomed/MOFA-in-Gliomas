source("setup.R")

## Load necessary data
legend_reactome_feature_set <- read.csv("https://reactome.org/download/current/ReactomePathways.txt", sep = "\t", header = FALSE)
legend_reactome_feature_set <- legend_reactome_feature_set[legend_reactome_feature_set$V3 == "Homo sapiens", ]
legend_reactome_feature_set <- legend_reactome_feature_set[,-3]
colnames(legend_reactome_feature_set) <- c("Pathway_id", "Pathway_name")

reactome_pathway_hierarchy <- read.csv("https://reactome.org/download/current/Complex_2_Pathway_human.txt", sep = "\t")



## Functions for GSEA 
run_rna_gsea <- function(loadings, sign, alpha){
  
  loadings <- sort(loadings, decreasing = T)
  ids_total <-  bitr(names(loadings), fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  ensembl_to_entrez <- setNames(ids_total$ENTREZID, ids_total$ENSEMBL)
  
  
  names(loadings) <- ensembl_to_entrez[names(loadings)]


  output <-  gsePathway(
    loadings,
    organism = "human",
    pvalueCutoff = alpha,
    pAdjustMethod = "BH",
    verbose = F,
    seed = T,
    by = "fgsea",
    minGSSize = 5,
    maxGSSize = 500
  )
  results <- output@result
  if (sign > 0){
    trunc_out <- results[results$NES > 0, ]
  }
  else if (sign < 0){
    trunc_out <- results[results$NES < 0, ]
  }
  else {
    trunc_out <- results
  }
  
  return(trunc_out)
}
  
  
run_methy_gsea <- function(factor, sign = c("all", "positive","negative"), background_vector, q, promoter = FALSE, alpha, filename = NULL){
  
  if (sign == "positive"){
    factor <- factor[factor > 0]
    q <- quantile(factor, 1-q)[[1]]
    factor_top <- factor[factor > q]
  }
  
  else if (sign == "negative") {
    factor <- factor[factor< 0]
    q <- quantile(factor, q)[[1]]
    factor_top <- factor[factor < q]
  }
  
  else {
    factor <- abs(factor)
    q <- quantile(factor, 1-q)[[1]]
    factor_top <- factor[factor > q]
  }
  sig_cgi <- names(factor_top)
  
  if (promoter){
    genomic.features_input = c("TSS1500","TSS200","1stExon")
  }
  else{
    genomic.features_input = "ALL"
  }
  
  output <- gometh(as.vector(sig_cgi),
                   all.cpg = background_vector,
                   collection =  "GO",  # c("GO", "KEGG"),
                   array.type = "450K", 
                   plot.bias = F,
                   prior.prob = TRUE,
                   anno = NULL,
                   equiv.cpg = TRUE,
                   fract.counts = TRUE,
                   genomic.features = genomic.features_input, 
                   sig.genes = FALSE
  ) 
  
  output <- output[output$ONTOLOGY == "BP", ]
  output <- output[order(output$P.DE, decreasing = FALSE),]
  output <- output[output$P.DE <= alpha & output$N >= 5 & output$N <=500, ]
  
  if (!is.null(filename)) 
    write.csv(as.data.frame(output), paste0("mofa-results/", filename, ".csv"), row.names = TRUE)
  
  return(output)
}


pathways_of_disease <- function(gsea_result){
  pathways <- gsea_result$ID
  
  sig_pathways_ids <- pathways
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
  paths_disease <- pathways[top_names == "Disease"]
  print(length(paths_disease))
  res <- c()
  for (id in paths_disease){
    v <- query(id = id)[["disease"]]
    res <- append(res, v$displayName)
  }
  return (res)
}




## Plot functions


plot_gsea <- function(gsea.results, alpha, max.pathways = 25, sign, filename = NULL) {
  
  if (sign > 0){
    gsea.results <- gsea.results[gsea.results$NES > 0, ]
    sign_color <-  "#268989"
  }
  else if (sign < 0){
    gsea.results <- gsea.results[gsea.results$NES < 0, ]
    sign_color <-  "#E43F3F"
  }
  else {
    gsea.results <- gsea.results
    sign_color <- "black"
  }
  
  
  gsea.results <- gsea.results[gsea.results$p.adjust <= alpha, ]

  # get p-values
  p.values <- gsea.results$p.adjust
  
  # Get data  
  tmp <- data.frame(
    pvalue = p.values, 
    pathway = gsea.results$ID
  )
  
  # if there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp)>max.pathways) tmp <- head(tmp[order(tmp$pvalue),], n=max.pathways)
  
  # convert pvalues to log scale
  tmp$logp <- -log10(tmp$pvalue+1e-100)
  
  # replace the pathways ids for the respective name
  tmp$pathway <- legend_reactome_feature_set$Pathway_name[match(tmp$pathway, legend_reactome_feature_set$Pathway_id)]
  tmp <- tmp[!is.na(tmp$pathway), ]
  rownames(tmp) <- tmp$pathway
  
  # order according to significance
  tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$start <- 0

 
  text_size <- 2
  p <- ggplot(tmp, aes(x=.data$pathway, y=.data$logp)) +
    geom_point(size=5,  color=sign_color) +
    geom_hline(yintercept=-log10(alpha), linetype="longdash") +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes(xend=.data$pathway, yend=.data$start)) +
    coord_flip() +
    labs(title = "", x = "", y = ("-log pvalue")) +
    my_theme
  
  if (!is.null(filename)) 
    ggsave(filename = paste0("mofa-results/", filename, ".pdf"), 
           plot = p, width = 14, height = 7, dpi = 300)
  
  
  return(p)
}


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
  
  p <- ggplot(plot_counts, aes(x = Factor, y = Count, fill = Top_Level_Pathway)) +
    geom_bar(stat = "identity", position = "fill") +  
    labs(title = "", 
         x = "", 
         y = "", 
         fill = "Top-Level Pathway") +
    scale_x_discrete(labels = names(list_sig_reactome_pathways_ids)) +
    my_theme
  
  if (!is.null(filename)) 
    ggsave(filename = paste0("mofa-results/", filename, ".pdf"), plot = p, width = 20, height = 7, dpi = 300)
  
  return(p)
}


plot_methy_gsea <- function(run_methy_gsea_output, sign, max.pathways = 25, filename = NULL,
                            alpha = 0.1){
  
  if (sign >0){
    sign_color <-  "#268989"
  }
  else if (sign <0) {
    sign_color <-  "#E43F3F"
  }
  else{
    sign_color <- "black"
  }
  
  up_table <- run_methy_gsea_output[,c("TERM","P.DE")][c(1: max.pathways), ]
  up_table$P.DE <- -log10(up_table$P.DE)
  up_table$TERM <- factor(up_table$TERM, levels = rev(up_table$TERM))
  
  p <- ggplot(up_table, aes(x=TERM, y=P.DE)) +
    geom_point(size=5, color=sign_color) +
    geom_hline(yintercept=-log10(alpha), linetype="longdash") +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes(xend=TERM, yend=0)) +
    my_theme +
    labs(x = "GO-BP term", y = expression(-log[10](p-value))) +
    coord_flip() 

  if (!is.null(filename)) 
    ggsave(filename = paste0("mofa-results/", filename, ".pdf"), plot = p, width = 14, height = 7, dpi = 300)
  
  return (p)
}



