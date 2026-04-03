setwd("~/Desktop/artigo")

source("setup.R")
### Load data
clinical <- read.csv("DataSets/2_analysis/clinical.csv", row.names = 1)

survival <- read.csv("DataSets/2_analysis/survival_data.csv", row.names = 1)

mutations <- read.csv("DataSets/2_analysis/mutations.csv", row.names = 1)
dna <- read.csv("DataSets/2_analysis/dna.csv", row.names = 1)
mrna <- read.csv("DataSets/2_analysis/mrna.csv", row.names = 1)
mirna <- read.csv("DataSets/2_analysis/mirna.csv", row.names = 1)


info.methy <- read.csv("DataSets/processed_assays/info_methylation.csv", row.names = 1)
info.mirna <- read.csv("DataSets/processed_assays/info_rna_mirna.csv", row.names = 1)
info.mrna <- read.csv("DataSets/processed_assays/info_rna_coding.csv", row.names = 1)

omics.list <- list(Mutations = mutations, Methylation = dna, mRNA = mrna, miRNA = mirna)




## Run Model
source("run_mofa.R")

output <- run_mofa(omics.list, n_runs = 10, num_factors = NULL, 
                     likelihoods = c("bernoulli", "gaussian", "gaussian", "gaussian"), 
                     clinical_data = clinical) 



#### Save Results
model <- output$model_output 
print(model)
print(output$time_taken)



#### Read Z and W 
factors <- as.matrix(read.csv("mofa-output/factors.csv", row.names = 1))

loadings <- list(
  Mutations = as.matrix(read.csv("mofa-output/loadings_mutations.csv", row.names = 1)),
  Methylation = as.matrix(read.csv("mofa-output/loadings_Methylation.csv", row.names = 1)),
  mRNA = as.matrix(read.csv("mofa-output/loadings_mrna.csv", row.names = 1)),
  miRNA = as.matrix(read.csv("mofa-output/loadings_mirna.csv", row.names = 1))
)

### Results shown in the paper 

## Variance Decomposition

p__variance_decomposition <- local({
  p1 <- plot_variance_explained(model, plot_total = T)[[1]]
  p1 <- p1$data
  
  p1$value <- round(p1$value, 1)
  colnames(p1) <- c("Factor", "View", "Value", "Group")
  p1$Factor <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", p1$Factor)
  p1$View <- factor(p1$View, levels = c("miRNA", "mRNA", "Methylation", "Mutations"))
  
  p <- ggplot(p1, aes(x = .data$Value, y = .data$View, fill = .data$Factor)) +
    geom_bar(
      position = position_stack(reverse = T), stat = "identity") +
    stat_identity(
      geom = "text", color = "black", size = 5,
      aes(label = ifelse(.data$Value > 2, .data$Value, "")),
      position = position_stack(reverse = TRUE, vjust = 0.5)
    ) +
    labs(title = "", x = "% of Variance", y = "") + 
    scale_fill_manual(values = brewer.pal(n = 4, name = "Oranges")) +
    my_theme
  
  p
})
ggsave(filename = "mofa-results/variance-decomposition.pdf",
       p__variance_decomposition,  width = 14, height = 7)



## Factors projection
classified_samples <- clinical[!is.na(clinical$Subtype), ]

p__factors_projection <- local({
  sample_df <- data.frame(Factor1 = numeric(nrow(classified_samples)),
                          Factor3 = numeric(nrow(classified_samples)),
                          Factor2 = numeric(nrow(classified_samples)),
                          Subtype = rep("GBM", nrow(classified_samples)))
  sample_df$Factor1 <- factors[rownames(classified_samples),1]
  sample_df$Factor3 <- factors[rownames(classified_samples),3]
  sample_df$Factor2 <- factors[rownames(classified_samples),2]
  sample_df$Subtype <- classified_samples$Subtype
  
  p <- ggplot(sample_df, aes(x = Factor1, y = Factor3)) +
    geom_point(aes(shape = Subtype, color = Factor2), size = 2.5) +
    geom_vline(xintercept=0, linetype="dashed") + 
    geom_segment(aes(x = 0, xend = Inf, y = 0, yend = 0), linetype = "dashed") +
    scale_color_gradientn(
      colours = c("darkgreen", "white", "darkred"),
      name = "Factor 2"
    ) +
    my_theme +
    labs(x = "Factor 1", y = "Factor 3", shape = "Subtype")
 p
})
ggsave(filename = "mofa-results/factors-projectionUPDATED.pdf", p__factors_projection, width = 6, height = 4)



## Hazard Ratio for factors
plot_survival_factors_univariate <- function(survival, factors_matrix) {
  
  survival[is.na(survival)] <- 0
  surv_object <- Surv(time = survival[rownames(factors_matrix), "Time"], 
                      event = survival[rownames(factors_matrix), "Status"])
  
  results <- lapply(colnames(factors_matrix), function(fac) {
    factor_values <- factors_matrix[, fac, drop = FALSE]
    fit <- coxph(surv_object ~ factor_values)
    s <- summary(fit)
    coef <- s$coefficients
    conf <- s$conf.int
    data.frame(
      factor = fac,
      coef = coef["factor_values", "exp(coef)"],
      p = coef["factor_values", "Pr(>|z|)"],
      lower = conf["factor_values", "lower .95"],
      higher = conf["factor_values", "upper .95"]
    )
  })
  
  df <- do.call(rbind, results)
  df$factor <- factor(df$factor, levels = rev(df$factor))
  df$significant <- ifelse(df$p < 0.05, "*", "")
  
  p <- ggplot(df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
    geom_pointrange() + 
    geom_text(aes(label = significant, y = coef), vjust = -1, size = 5)+
    labs(y="Hazard Ratio", x="") + 
    scale_x_discrete(labels = paste0("Factor ", seq(ncol(factors_matrix),1))) + 
    geom_hline(aes(yintercept=1), linetype="dotted") +
    coord_flip() +
    my_theme
  
  return(list(p, df))
}

plot_survival_factors_multivariate <- function (survival, factors_matrix){
  
  survival[is.na(survival)] <- 0
  surv_object <- Surv(time = survival[rownames(factors_matrix), "Time"], event = survival[rownames(factors_matrix), "Status"] )
  
  fit <- coxph(surv_object ~ factors_matrix) 
  
  s <- summary(fit)
  coef <- s[["coefficients"]]
  
  df <- data.frame(
    factor = factor(rownames(coef), levels = rev(rownames(coef))),
    p      = coef[,"Pr(>|z|)"], 
    coef   = coef[,"exp(coef)"], 
    lower  = s[["conf.int"]][,"lower .95"], 
    higher = s[["conf.int"]][,"upper .95"])
  
  df$significant <- ifelse(df$p < 0.05, "*", "")
  
  p <- ggplot(df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
    geom_pointrange() + 
    geom_text(aes(label = significant, y = coef), vjust = -1, size = 5)+
    labs(y="Hazard Ratio", x="") + 
    scale_x_discrete(labels = paste0("Factor ", seq(ncol(factors_matrix),1))) + 
    geom_hline(aes(yintercept=1), linetype="dotted") +
    coord_flip() +
    my_theme
  print(p)
  
  return(list(p, fit))
}

p <- plot_survival_factors_univariate(survival, factors)[[1]]
ggsave(filename = "mofa-results/factors-hazard-uni.pdf", p, width = 6, height = 4)

p <- plot_survival_factors_multivariate(survival, factors)[[1]]
ggsave(filename = "mofa-results/factors-hazard-multi.pdf", p, width = 6, height = 4)

rm(p)


## Features selected & mutations projections

# Mutations Omic
mutations_sel_f1_pos <- c("IDH1")
mutations_sel_f1_neg <- c("PTEN", "EGFR")

mutations_sel_f3_pos <- c("TP53", "ATRX")
mutations_sel_f3_neg <- c("CIC")

p__mutations_projection <- local({
  combs <- c("EGFR", "IDH1", "IDH1", "ATRX", "CIC", "TP53")
  p_list <- list()
  for (k in c(1, 3, 5)){
    p <- plot_factors(model, 
                      factors = c(1,3), 
                      color_by = combs[k],
                      shape_by = combs[k+1],
                      show_missing = T,
                      scale = T
    )
    p <- p + 
      geom_vline(xintercept=0, linetype="dashed") + 
      geom_segment(aes(x = 0, xend = max(p$data$x), y = 0, yend = 0), linetype = "dashed") +
      my_theme +
      labs(x = "Factor 1", y = "Factor 3")
    
    ggsave(filename = paste0("mofa-results/mutations-projected-", k, ".pdf"), p, width = 6, height = 4)
    p_list <- append(p_list, p)
  }
  
  p_list
})
  

# DNA Methylation Omic
methy_sel_f1 <- names(head(sort(abs(loadings$Methylation[, 1]), decreasing = TRUE), 30)) 
methy_sel_f3 <- names(head(sort(abs(loadings$Methylation[, 3]), decreasing = TRUE), 30))

# mRNA Omic
mrna_sel_f1 <- names(head(sort(abs(loadings$mRNA[, 1]), decreasing = TRUE), 30))
mrna_sel_f2 <- names(head(sort(abs(loadings$mRNA[, 2]), decreasing = TRUE), 30))
mrna_sel_f3 <- names(head(sort(abs(loadings$mRNA[, 3]), decreasing = TRUE), 30))

# miRNA Omic
mirna_sel_f1 <- names(head(sort(abs(loadings$miRNA[, 1]), decreasing = TRUE), 3))



## Enrichment Analysis
source("enrichment.R")

# mRNA
mrna_enrich_f1 <- run_rna_gsea(loadings$mRNA[, 1], sign = 0, alpha = 0.05)
plot_gsea(mrna_enrich_f1, alpha = 0.05, max.pathways = 25, sign = -2, filename = NULL)

mrna_enrich_f2 <- run_rna_gsea(loadings$mRNA[, 2], sign = 0, alpha = 0.05)
plot_gsea(mrna_enrich_f2, alpha = 0.05, max.pathways = 25, sign = 2, filename = NULL)

mrna_enrich_f3 <- run_rna_gsea(loadings$mRNA[, 3], sign = 0, alpha = 0.05)
plot_gsea(mrna_enrich_f3, alpha = 0.05, max.pathways = 25, sign = 2, filename = NULL)

plot_categories_reactome(list("Factor 1"= mrna_enrich_f1[mrna_enrich_f1$NES<0,"ID"], 
                              "Factor 2"= mrna_enrich_f2[mrna_enrich_f2$NES>0, "ID"], 
                              "Factor 3"= mrna_enrich_f3[mrna_enrich_f3$NES>0, "ID"]), 
                         filename="mrna-gsea-top-levels") 
#the pathways in "Disease"
diseases_f1 <- pathways_of_disease(mrna_enrich_f1)
print(table(diseases_f1))

diseases_f2 <- pathways_of_disease(mrna_enrich_f2)
print(table(diseases_f2))

diseases_f3 <- pathways_of_disease(mrna_enrich_f3)
print(table(diseases_f3))


# DNA methylation
methy_enrich_f1 <- run_methy_gsea(loadings$Methylation[,1], sign = "positive", 
                                  background_vector = rownames(loadings$Methylation), 
                                  q = 0.01, promoter = F, alpha = 0.05, filename = NULL)
plot_methy_gsea(methy_enrich_f1, sign = 2, max.pathways = 15, 
                      filename = "methy-enrich-f1", 
                      alpha = 0.05)

methy_enrich_f1_promoter <- run_methy_gsea(loadings$Methylation[,1], sign = "positive", 
                                           background_vector = rownames(loadings$Methylation), q = 0.01, 
                                           promoter = TRUE, alpha = 0.05, filename = NULL)

plot_methy_gsea(methy_enrich_f1_promoter, sign = 2, max.pathways = 15, 
                      filename = "methy-enrich-f1-prom", 
                      alpha = 0.05)



methy_enrich_f3 <- run_methy_gsea(loadings$Methylation[,3], sign = "negative", 
                                  background_vector = rownames(loadings$Methylation), 
                                  q = 0.01, promoter = F, alpha = 0.05, filename = NULL)
plot_methy_gsea(methy_enrich_f3, sign = -2, max.pathways = 15, 
                filename = "methy-enrich-f3", 
                alpha = 0.05)


methy_enrich_f3_promoter <- run_methy_gsea(loadings$Methylation[,3], sign = "negative", 
                                           background_vector = rownames(loadings$Methylation), 
                                           q = 0.01, 
                                           promoter = TRUE, alpha = 0.05, filename = NULL)
plot_methy_gsea(methy_enrich_f3_promoter, sign = -2, max.pathways = 15, 
                filename = "methy-enrich-f3-prom", 
                alpha = 0.05)





## Validation Features


# Survival
gbm_mrna <- read.csv("survival-results/gbm_mrna.csv", row.names = 1)
lgg_mrna <- read.csv("survival-results/lgg_mrna.csv", row.names = 1)

gbm_methy <- read.csv("survival-results/gbm_methy.csv", row.names = 1)
lgg_methy <- read.csv("survival-results/lgg_methy.csv", row.names = 1)



# DGE
gbm_vs_astro_mrna <- read.csv("dge-results/gbm_vs_astro_mrna.csv", row.names = 1)
gbm_vs_oligo_mrna <- read.csv("dge-results/gbm_vs_oligo_mrna.csv", row.names = 1)
astro_vs_oligo_mrna <- read.csv("dge-results/astro_vs_oligo_mrna.csv", row.names = 1)

gbm_vs_astro_methy <- read.csv("dge-results/gbm_vs_astro_methy.csv", row.names = 1)
gbm_vs_oligo_methy <- read.csv("dge-results/oligo_vs_gbm_methy.csv", row.names = 1)
astro_vs_oligo_methy <- read.csv("dge-results/astro_vs_oligo_methy.csv", row.names = 1)

# build data frame as summary for each factor:
                              #DGE (GBM vs ASTRO; GBM vs OLIGO; ASTRO vs OLIGO);
                              #Significant in survival: in GBM group || in LGG group
build_df <- function(methy_sel, mrna_sel) {
  
  combined_features <- c(methy_sel, mrna_sel)
  df <- data.frame(matrix(NA, nrow = length(combined_features), ncol = 5))
  colnames(df) <- c("GBMvsASTRO", "GBMvsOLIGO", "ASTROvsOLIGO",  "GBM", "LGG")
  rownames(df) <- combined_features
  
  df$"GBMvsASTRO" <- c(gbm_vs_astro_methy[methy_sel, "diffexpressed"],
                       gbm_vs_astro_mrna[mrna_sel, "diffexpressed"])
  
  df$"GBMvsOLIGO" <- c(gbm_vs_oligo_methy[methy_sel, "diffexpressed"],
                       gbm_vs_oligo_mrna[mrna_sel, "diffexpressed"])
  
  df$"ASTROvsOLIGO" <- c(astro_vs_oligo_methy[methy_sel, "diffexpressed"],
                         astro_vs_oligo_mrna[mrna_sel, "diffexpressed"])
  
  df$"LGG" <- round(c(lgg_methy[methy_sel, "P_Value"],
                lgg_mrna[mrna_sel, "P_Value"]),3)
  
  df$"GBM" <- round(c(gbm_methy[methy_sel, "P_Value"],
                gbm_mrna[mrna_sel, "P_Value"]),3)
  return(df)
}
df_f1 <- build_df(methy_sel_f1, mrna_sel_f1)
df_f2 <- build_df(character(0), mrna_sel_f2)
df_f3 <- build_df(methy_sel_f3, mrna_sel_f3)

write.csv(df_f1,"summary-f1.csv", row.names = TRUE)
write.csv(df_f2,"summary-f2.csv", row.names = TRUE)
write.csv(df_f3,"summary-f3.csv", row.names = TRUE)



## Correlation of CpG and mRNA
plot_cor_cpgs_mrna <- local({
  conc_matrix <- cbind(omics.list$Methylation, omics.list$mRNA)
  conc_matrix <- conc_matrix[,c(methy_sel_f1, mrna_sel_f1, methy_sel_f3, mrna_sel_f3)]
  correlation_matrix <- cor(conc_matrix, method = "pearson", use = "na.or.complete")
  rownames(correlation_matrix) <- c(methy_sel_f1, 
                                    info.mrna[rownames(info.mrna) %in% mrna_sel_f1, "gene_name"], 
                                    methy_sel_f3, 
                                    info.mrna[rownames(info.mrna) %in%mrna_sel_f3, "gene_name"])
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#268989", "white", "#E43F3F"))
  
  ht <- Heatmap(
    correlation_matrix,
    name = "Pearson Correlation",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 20),
    heatmap_legend_param = list(
      title = NULL,  
      title_gp = gpar(fontsize = 0),  
      labels_gp = gpar(fontsize = 20),
      legend_height = unit(75, "cm")  
    )
  )
  pdf("mofa-results/heatmap-cpgs-mrna.pdf", width = 30, height = 30)
  draw(ht)
  dev.off()
})



## Cpg-Gene regulation
#in particular, gene ISM1 with the 8 probes selected
cpgs_f3_ism1 <- intersect(methy_sel_f3, rownames(info.methy)[info.methy$gene=="ISM1"])
id_ism1 <- rownames(info.mrna)[info.mrna$gene_name=="ISM1"]
expr_ism1 <- mrna[, id_ism1, drop = FALSE]

lgg_patients <- rownames(clinical)[clinical$Type == "LGG" & !is.na(clinical$Type)]
gbm_patients <- rownames(clinical)[clinical$Type == "GBM" & !is.na(clinical$Type)]

# correlations of each cpg with the gene ISM1
cors <- local({
  cors <- c()
  for (cpg in cpgs_f3_ism1){
    val <- cor(expr_ism1[lgg_patients, ], dna[lgg_patients, cpg], method = "pearson", use = "complete.obs")
    cors <- append(cors, val)
  }
  cors
})
cors

# correlations of each cpg with the genes selected by factor 3 (choose idx)
cors_matrix <- local({
  cors_matrix <- matrix(data = NA, 
                        nrow = length(cpgs_f3_ism1),
                        ncol = length(mrna_sel_f3),
                        dimnames = list(cpgs_f3_ism1, mrna_sel_f3))
  for (idx_cpg in 1: length(cpgs_f3_ism1)){
    for (idx_gene in 1: length(mrna_sel_f3)){
      gene_name <- info.mrna[rownames(info.mrna) == mrna_sel_f3[idx_gene], "gene_name"]
      val <- cor(mrna[lgg_patients, mrna_sel_f3[idx_gene]], dna[lgg_patients, cpgs_f3_ism1[idx_cpg]], 
                 method = "pearson", use = "complete.obs")
      cors_matrix[idx_cpg, idx_gene] <- val 
    }
  }
  colnames(cors_matrix) <- info.mrna[rownames(info.mrna) %in% mrna_sel_f3, "gene_name"]
  rownames(cors_matrix) <- cpgs_f3_ism1
  cors_matrix
})


# Plot
cor_plot_gene_cpg <- function(cpg, gene_id, labels, groups_2_test){
  idxs <- which(labels %in% groups_2_test) 
  expr_gene <- mrna[idxs , gene_id, drop = FALSE]
  gene_name <- info.mrna[rownames(info.mrna) == gene_id, "gene_name"]
  expr_cpg <- dna[idxs, cpg, drop = FALSE]
  
  df_plot <- data.frame(
    expr_gene = expr_gene[,1],
    expr_cpg = expr_cpg[,1],
    groups = labels[idxs]
  )
  
  plot <- ggplot(df_plot, aes(x = expr_cpg, y = expr_gene, color = groups)) +
    geom_smooth(method = "lm", se = T) +
    geom_point() +
    scale_color_manual(values = c('#FF7400', '#009999', "black")) +
    stat_cor(method = "pearson",
             cor.coef.name = "R",
             label.x = 1.08) +
    labs(
      x = cpg,,
      y =  gene_name,
      color = "Subtype"
    ) + 
    my_theme
  ggsave(filename = paste0("mofa-results/corr-",cpg, "-", gene_name,".pdf"), plot, width = 8, height = 6)
  return (plot)
}

p <- cor_plot_gene_cpg(cpgs_f3_ism1[1], id_ism1, 
                       labels = clinical$Subtype, groups_2_test = c("ASTRO", "OLIGO", "GBM"))
p <- cor_plot_gene_cpg(cpgs_f3_ism1[2], id_ism1, 
                       labels = clinical$Subtype, groups_2_test = c("ASTRO", "OLIGO", "GBM"))

# heatmap of cors_matrix
col_fun <- colorRamp2(c(-1, 0, 1), c("#268989", "white", "#E43F3F"))

cors_matrix <- cbind(cors_matrix, "ISM1" = cors)

ht <- Heatmap(
  cors_matrix,
  name = "Pearson Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
    show_column_names = T,
  row_names_gp = gpar(fontsize = 20),
  heatmap_legend_param = list(
    title = NULL,  
    title_gp = gpar(fontsize = 0),  
    labels_gp = gpar(fontsize = 20),
    legend_height = unit(35, "cm")  
  )
)
pdf("mofa-results/heatmap-mrnaf3-cpgsiNism1.pdf", width = 20, height = 15)
draw(ht)
dev.off()


#  SLC2A5       BLNK  TMEM119 PLXDC2
id_plxdc2 <- rownames(info.mrna)[info.mrna$gene_name=="PLXDC2"]
id_slc2a5 <- rownames(info.mrna)[info.mrna$gene_name=="SLC2A5"]
id_blnk <- rownames(info.mrna)[info.mrna$gene_name=="BLNK"]
id_tmem119 <- rownames(info.mrna)[info.mrna$gene_name=="TMEM119"]

p <- cor_plot_gene_cpg(cpgs_f3_ism1[1], id_plxdc2, 
                       labels = clinical$Subtype, groups_2_test = c("ASTRO", "OLIGO", "GBM"))
p <- cor_plot_gene_cpg(cpgs_f3_ism1[2], id_plxdc2, 
                       labels = clinical$Subtype, groups_2_test = c("ASTRO", "OLIGO", "GBM"))


