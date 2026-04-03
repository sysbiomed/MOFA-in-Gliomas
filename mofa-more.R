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



#### Read Z and W 
factors <- as.matrix(read.csv("mofa-output/factors.csv", row.names = 1))

loadings <- list(
  Mutations = as.matrix(read.csv("mofa-output/loadings_mutations.csv", row.names = 1)),
  Methylation = as.matrix(read.csv("mofa-output/loadings_Methylation.csv", row.names = 1)),
  mRNA = as.matrix(read.csv("mofa-output/loadings_mrna.csv", row.names = 1)),
  miRNA = as.matrix(read.csv("mofa-output/loadings_mirna.csv", row.names = 1))
)



### Find the best number of clusters
distribution_logp <- local({
  res_l <- list()
  clust_l <- list()
  for (k in 3: 8) {  
    set.seed(4)  # for reproducibility
    res <- c()
    
    for (r in 1: 30) {
      strat_sample <- survival %>%
        tibble::rownames_to_column("SampleID") %>%
        group_by(Status) %>%
        sample_frac(size = 0.80) %>%
        tibble::column_to_rownames("SampleID")
      
      km_res <- kmeans(factors[rownames(strat_sample), c(1,2,3)], centers = k, nstart = 10)
  
      surv_obj <- Surv(time = strat_sample[, "Time"] , 
                       event = strat_sample[, "Status"])
      
      fit <- survfit(surv_obj ~ factor(km_res$cluster))
      logrank <- survdiff(surv_obj ~ factor(km_res$cluster))
      
      p_val <- logrank$pvalue
      res <- c(res, p_val)
    }
    res_l[[as.character(k)]] <- res
  }
  res_l
})

plot_distribution_logp <- local({
  res_l_log10 <- sapply(distribution_logp, function(x) -log10(x))
  
  df <- data.frame(res_l_log10)
  colnames(df) <- c(3,4,5,6,7,8)
  
  df_long <- df %>%
    pivot_longer(
      cols = everything()
    )
  
  p <- ggplot(df_long, aes(x = name, y = value, fill = name)) +
    geom_boxplot(fill = "#FF7F32") +
    labs(x = "Number of Clusters", y = expression(-log[10](p-value))) +
    my_theme
  
  ggsave(filename = paste0("mofa-more/", "distribution-logp-k", ".pdf"), 
         plot = p, width = 6, height = 4, dpi = 300)  
})


### Functions

plot_projection <- function(km_res){
  k <- nrow(km_res$centers)
  plot_proj_data <- data.frame(x = factors[,1], 
                               y  = factors[,3], 
                               color_by = factor(km_res$cluster,
                                                 levels = 1:k, 
                                                 labels = paste("Cluster", 1:k)))
  
  
  my_colors <- brewer.pal(n = k, name = "Dark2")
  p <- ggplot(plot_proj_data , aes(x=.data$x, y=.data$y, color=.data$color_by)) + 
    geom_point(size=3.0) +
    labs(x="Factor 1", y="Factor 3", color = "") +
    scale_color_manual(values = my_colors) +
    my_theme +
    guides(color = guide_legend(title = "Cluster"))
  
  ggsave(filename = paste0("mofa-more/projection-", k , ".pdf"), plot = p, width = 10, height = 7, dpi = 300)
  
}


set.seed(44)

### Results for k = 3

km_res_3 <- kmeans(factors[,c(1,2,3)], centers = 3, nstart = 10)

# Projection

plot_projection(km_res_3)

# Survival Curves

survival_curves_3 <- local({
  plot_data <- data.frame(
    Time = survival[, "Time"],
    Status = survival[ ,"Status"],
    Cluster = factor(km_res_3$cluster,levels = 1:3, 
                     labels = paste("Cluster", 1:3))
  )
  surv_obj <- Surv(time = plot_data$Time, event = plot_data$Status)
  fit <- survfit(surv_obj ~ Cluster, data = plot_data)
  names(fit$strata) <- gsub("Cluster=", "", names(fit$strata))
  
  pairwise_res <- pairwise_survdiff(Surv(Time, Status) ~ Cluster, data = plot_data, p.adjust.method = "BH")
  pairwise_pvals <- as.data.frame(pairwise_res$p.value)
  pval_text <- paste0("Clusters 1 vs 2: ", signif(pairwise_pvals[1,1], 3),"\n",
                      "Clusters 1 vs 3: ", signif(pairwise_pvals[2,1], 3), "\n",
                      "Clusters 2 vs 3: ", signif(pairwise_pvals[2,2], 3), "\n")
  
  p <- ggsurvplot(
    fit,
    data = plot_data,
    risk.table = FALSE,
    pval = FALSE,
    palette = "Dark2",
    legend.title = "",
    xlab = "Time (days)",
    ylab = "Survival Probability"
  )
  
  # Add annotation with p-values
  p$plot <- p$plot +
    annotate("text", x = 3000, y = 0.8, label = pval_text, hjust = 0, size = 6) +
    my_theme
  p
  ggsave(filename = "mofa-more/surv-curves-3.pdf", plot = p$plot, width = 10, height = 8, dpi = 300)
})




# Table of sample information
summary_cluster_3 <- local({
  table_3 <- clinical
  table_3$Cluster <- factor(km_res_3$cluster,levels = 1:3, 
                            labels = paste("Cluster", 1:3))
  table_3$IDH1 <- omics.list$Mutations[, "IDH1"]
  table_3$EGFR <- omics.list$Mutations[, "EGFR"]
  table_3$PTEN <- omics.list$Mutations[, "PTEN"]
  table_3$TP53 <- omics.list$Mutations[, "TP53"]
  table_3$PIK3CA <- omics.list$Mutations[, "PIK3CA"]
  
  cluster_summary <- table_3 %>%
    group_by(Cluster) %>%
    summarise(
      n = n(),
      AGE = round(mean(Age, na.rm = TRUE), 2),
      MALE = round(mean(Sex == "Male", na.rm = TRUE) * 100,2),
      OLIGO = sum(Subtype == "OLIGO", na.rm = TRUE),
      ASTRO = sum(Subtype == "ASTRO", na.rm = TRUE),
      GBM = sum(Subtype == "GBM", na.rm = TRUE),
      DEAD = round(sum(Vital_Status == "Dead", na.rm = TRUE)/n * 100,2),
      IDH1 = round(sum(IDH1 == 1, na.rm = TRUE)/n * 100,2),
      PTEN = round(sum(PTEN == 1, na.rm = TRUE)/n * 100,2),
      EGFR = round(sum(EGFR == 1, na.rm = TRUE)/n * 100,2),
      TP53 = round(sum(TP53 == 1, na.rm = TRUE)/n * 100,2),
      PIK3CA = round(sum(PIK3CA == 1, na.rm = TRUE)/n * 100,2)
    )
  cluster_summary
  write.csv(cluster_summary, "mofa-more/sample-info-3.csv", row.names = TRUE)
})



### Results for k = 5
km_res_5 <- kmeans(factors[,c(1,2,3)], centers = 5, nstart = 10)

# Projection
plot_projection(km_res_5)

# Survival Curves
survival_curves_5 <- local({
  
  plot_data <- data.frame(
    Time = survival[, "Time"],
    Status = survival[ ,"Status"],
    Cluster = factor(km_res_5$cluster,levels = 1:5, 
                     labels = paste("Cluster", 1:5))
  )
  
  surv_obj <- Surv(time = plot_data[,"Time"], event = plot_data[,"Status"])
  fit <- survfit(surv_obj ~ Cluster, data = plot_data)
  names(fit$strata) <- gsub("Cluster=", "", names(fit$strata))
  
  p <- ggsurvplot(
    fit,
    data = plot_data,
    risk.table = F,
    pval = F,
    palette = "Dark2",
    legend.title = "",
    xlab = "Time (days)",
    ylab = "Survival Probability"
  )
  
  p <- p$plot + 
    my_theme 
  
  ggsave(filename = "mofa-more/surv-curves-5.pdf", plot = p, width = 10, height = 8, dpi = 300)

  pairwise_res <- pairwise_survdiff(Surv(Time, Status) ~ Cluster, data = plot_data, p.adjust.method = "BH")
  
  # Extract p-values into a readable format
  pairwise_pvals <- as.data.frame(pairwise_res$p.value)
  pairwise_pvals_rounded <- signif(pairwise_pvals, digits=3)
  print(pairwise_pvals_rounded)
  write.csv(pairwise_pvals_rounded, "mofa-more/pairwise-pval-5.csv", row.names = TRUE)
})



# Table of sample information
summary_cluster_5 <- local({
  table_5 <- clinical
  table_5$Cluster <- factor(km_res_5$cluster,levels = 1:5, 
                            labels = paste("Cluster", 1:5))
  table_5$IDH1 <- omics.list$Mutations[, "IDH1"]
  table_5$EGFR <- omics.list$Mutations[, "EGFR"]
  table_5$PTEN <- omics.list$Mutations[, "PTEN"]
  table_5$TP53 <- omics.list$Mutations[, "TP53"]
  table_5$PIK3CA <- omics.list$Mutations[, "PIK3CA"]
  
  cluster_summary <- table_5 %>%
    group_by(Cluster) %>%
    summarise(
      n = n(),
      AGE = round(mean(Age, na.rm = TRUE), 2),
      MALE = round(mean(Sex == "Male", na.rm = TRUE) * 100,2),
      OLIGO = sum(Subtype == "OLIGO", na.rm = TRUE),
      ASTRO = sum(Subtype == "ASTRO", na.rm = TRUE),
      GBM = sum(Subtype == "GBM", na.rm = TRUE),
      DEAD = round(sum(Vital_Status == "Dead", na.rm = TRUE)/n * 100,2),
      IDH1 = round(sum(IDH1 == 1, na.rm = TRUE)/n * 100,2),
      PTEN = round(sum(PTEN == 1, na.rm = TRUE)/n * 100,2),
      EGFR = round(sum(EGFR == 1, na.rm = TRUE)/n * 100,2),
      TP53 = round(sum(TP53 == 1, na.rm = TRUE)/n * 100,2),
      PIK3CA = round(sum(PIK3CA == 1, na.rm = TRUE)/n * 100,2)
    )
  cluster_summary
  write.csv(cluster_summary, "mofa-more/sample-info-5.csv", row.names = TRUE)
})


# Histological Type
classification <- as.matrix(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", row.names = 1))
missing_samples <- setdiff(rownames(omics.list$Mutations),rownames(classification))
new_entries <- matrix("unclassified", 
                      nrow = length(missing_samples), 
                      ncol = ncol(classification),
                      dimnames = list(missing_samples, colnames(classification)))


classification <- rbind(classification, new_entries)


classification_1 <- classification[names(km_res_5$cluster)[km_res_5$cluster == 1], ]
classification_2 <- classification[names(km_res_5$cluster)[km_res_5$cluster == 2], ]
classification_3 <- classification[names(km_res_5$cluster)[km_res_5$cluster == 3], ]
classification_4 <- classification[names(km_res_5$cluster)[km_res_5$cluster == 4], ]
classification_5 <- classification[names(km_res_5$cluster)[km_res_5$cluster == 5], ]

table(classification_1[, "TCGA.histological.type"]) 
table(classification_5[, "TCGA.histological.type"]) 
table(classification_3[, "TCGA.histological.type"])
table(classification_2[, "TCGA.histological.type"]) 
table(classification_4[, "TCGA.histological.type"])


### Re-label
clust_labeled <- as.matrix(km_res_5$cluster)


new_labels <- c(
  "1" = "GBM_1",
  "2" = "ASTRO_2",
  "3" = "GBM_3",
  "4" = "MIX_LGG_4",
  "5" = "OLIGO_5"
)
clust_labeled_named <- new_labels[as.character(clust_labeled[, 1])]
names(clust_labeled_named) <- rownames(clust_labeled)
clust_labeled_named <- as.vector(clust_labeled_named)


### Differences between clusters (DGE and DMP) - for k = 5


## DGE
mrna <- as.matrix(read.csv("DataSets/original_assays/assay_rna_coding.csv", row.names = 1))

mrna <- t(mrna)

keep <- rowSums(cpm(mrna) > 1) >= 5
mrna <- mrna[keep,]

mrna <- DGEList(counts = mrna)
mrna <- normLibSizes(mrna, method = "TMM")

patients <- rownames(omics.list$Mutations)

# only keep the 318 patients
mrna$samples <- mrna$samples[patients, ]
mrna$counts <- mrna$counts[,patients]

mrna$samples$group <- clust_labeled_named

group <- factor(clust_labeled_named)

design <- model.matrix(~ 0 + group, data = mrna$samples)
colnames(design) <- levels(group)
mrna <- estimateDisp(mrna, design, robust=TRUE)
fit <- glmQLFit(mrna, design, robust=TRUE)


group_names <- c("GBM_1", "ASTRO_2", "GBM_3", "MIX_LGG_4", "OLIGO_5")
contrast_pairs <- combn(group_names, 2, simplify = FALSE)

contrast_results <- list()

for (pair in contrast_pairs) {
  contrast_name <- paste0(pair[1], "_vs_", pair[2])
  contrast_expr <- paste0(pair[1], " - ", pair[2])
  
  contrast <- makeContrasts(contrasts = contrast_expr, levels = design)
  fit_result <- glmQLFTest(fit, contrast = contrast)
  
  contrast_results[[contrast_name]] <- fit_result
}

dge_2_save_mrna <- function(df){
  
  df$diffexpressed <- "NO"
  df$diffexpressed[!is.na(df$logFC)  & df$logFC > 1 & df$PValue < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$PValue < 0.05] <- "DOWN"
  
  return (df)
}


dfs_mrna <- lapply(contrast_results, function(x) dge_2_save_mrna(x$table))
lapply(names(dfs_mrna), function(name) {
  write.csv(dfs_mrna[[name]],
            file = paste0("mofa-more/dge-", name, ".csv"))
})



## DMP
dna <- as.matrix(read.csv("DataSets/processed_assays/assay_methylation.csv", row.names = 1))
dna <-  dna[rownames(omics.list$Mutations), ]


dna_beta <- t(M2B(dna))

unique(clust_labeled_named)


#factor 1
fac_1_de_dna <- matrix(data = NA, nrow = 5, ncol = 5, dimnames = list(new_labels, new_labels))
for (i in 1:5){
  for (j in (i+1):5){
    contrast_str <- c(colnames(fac_1_de_dna)[j], colnames(fac_1_de_dna)[i])
    
    res <- champ.DMP(beta = dna_beta,
                     pheno = clust_labeled_named,
                     compare.group = contrast_str,
                     adjPVal = 0.05 ,
                     adjust.method = "BH",
                     arraytype = "450K")
    res <- res[[1]]
    res <- res[methy_sel_f1, ]
    
    de_probes <- rownames(res)[abs(res$logFC) > 0.3 &
                                 res$adj.P.Val < 0.05
    ]
    
    fac_1_de_dna[i, j] <- length(de_probes)
  }
}
fac_1_de_dna[lower.tri(fac_1_de_dna)] <- t(fac_1_de_dna)[lower.tri(fac_1_de_dna)]
print(fac_1_de_dna)


# factor 3
fac_3_de_dna <- matrix(data = NA, nrow = 5, ncol = 5, dimnames = list(new_labels, new_labels))

for (i in 1:5){
  for (j in (i+1):5){
    contrast_str <- c(colnames(fac_3_de_dna)[j], colnames(fac_3_de_dna)[i])
    
    res <- champ.DMP(beta = dna_beta,
                     pheno = clust_labeled_named,
                     compare.group = contrast_str,
                     adjPVal = 0.05 ,
                     adjust.method = "BH",
                     arraytype = "450K")
    res <- res[[1]]
    res <- res[methy_sel_f3, ]
    
    de_probes <- rownames(res)[abs(res$logFC) > 0.3 &
                                 res$adj.P.Val <= 0.05
    ]
    
    fac_3_de_dna[i, j] <- length(de_probes)
  }
}
fac_3_de_dna[lower.tri(fac_3_de_dna)] <- t(fac_3_de_dna)[lower.tri(fac_3_de_dna)]
print(fac_3_de_dna)