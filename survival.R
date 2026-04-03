source("setup.R")

## Load data 
mrna <- as.matrix(read.csv("DataSets/processed_assays/assay_rna_coding.csv", row.names = 1))
dna <- as.matrix(read.csv("DataSets/processed_assays/assay_methylation.csv", row.names = 1))

classification <- as.matrix(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", row.names = 1))
survival.complete <- as.matrix(read.csv("DataSets/processed_assays/survival_complete.csv", row.names = 1))


classification <- classification[,"classification.2021"]
classification <- ifelse(classification == "glioblastoma", "GBM", 
                         ifelse(classification  == "astrocytoma", "ASTRO",
                                ifelse(classification  == "oligodendroglioma", "OLIGO",
                                       NA)))

loadings <- list(
  Mutations = as.matrix(read.csv("mofa-output/loadings_mutations.csv", row.names = 1)),
  Methylation = as.matrix(read.csv("mofa-output/loadings_Methylation.csv", row.names = 1)),
  mRNA = as.matrix(read.csv("mofa-output/loadings_mrna.csv", row.names = 1)),
  miRNA = as.matrix(read.csv("mofa-output/loadings_mirna.csv", row.names = 1))
)

# DNA Methylation Omic
methy_sel_f1 <- names(head(sort(abs(loadings$Methylation[, 1]), decreasing = TRUE), 30)) 
methy_sel_f3 <- names(head(sort(abs(loadings$Methylation[, 3]), decreasing = TRUE), 30))

# mRNA Omic
mrna_sel_f1 <- names(head(sort(abs(loadings$mRNA[, 1]), decreasing = TRUE), 30))
mrna_sel_f2 <- names(head(sort(abs(loadings$mRNA[, 2]), decreasing = TRUE), 30))
mrna_sel_f3 <- names(head(sort(abs(loadings$mRNA[, 3]), decreasing = TRUE), 30))


## Arrange data
mrna <- mrna[intersect(rownames(mrna), intersect(rownames(survival.complete), names(classification))), ] # 652
dna <- dna[intersect(rownames(dna), intersect(rownames(survival.complete), names(classification))), ] # 624

survival.complete <- survival.complete[rownames(mrna), ]
classification <- classification[rownames(mrna)]
survival.complete <- as.data.frame(survival.complete)

survival.complete_dna <- survival.complete[rownames(dna), ]
classification_dna <- classification[rownames(dna)]
survival.complete_dna <- as.data.frame(survival.complete_dna)


## Function

surv_group <- function(vec_genes, omic_dataset, group, survival, classes){
  
  res <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  if (length(group) > 1){
    samples <- c(names(classes == group[1]), names(classes == group[2]))
  }
  else{
    samples <- c(names(classes == group))
  }
  data <- survival[samples, ]
  
  for (gene in vec_genes) {
    median <- median(omic_dataset[samples, gene], na.rm = TRUE)
    data$gene_group <- ifelse(omic_dataset[samples, gene] >= median, "High", "Low")
    
    surv_obj <- Surv(time = survival.complete[samples, "Time"], 
                     event = survival.complete[samples, "Status"])
    log_rank_test <- survdiff(surv_obj ~ gene_group, data = data)
    p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1) 
    
    res <- rbind(res, data.frame(Gene = gene, P_Value = p_value))
  }
  rownames(res) <- res[,1]
  res <- res[, -1, drop = F]
  return (res)
}



gbm_mrna <- surv_group(c(mrna_sel_f1, mrna_sel_f2, mrna_sel_f3), mrna, "GBM",
                       survival.complete, classification)
lgg_mrna <- surv_group(c(mrna_sel_f1, mrna_sel_f2, mrna_sel_f3), mrna, c("ASTRO", "OLIGO"), 
                       survival.complete, classification)
                       
gbm_methy <- surv_group(c(methy_sel_f1, methy_sel_f3), dna, "GBM", 
                              survival.complete_dna, classification_dna)
lgg_methy <- surv_group(c(methy_sel_f1, methy_sel_f3), dna, c("ASTRO", "OLIGO"),
                              survival.complete_dna, classification_dna)
                                              
                       
                       

## SAVE
write.csv(gbm_mrna, "survival-results/gbm_mrna.csv", row.names = TRUE)
write.csv(lgg_mrna, "survival-results/lgg_mrna.csv", row.names = TRUE)

write.csv(gbm_methy, "survival-results/gbm_methy.csv", row.names = TRUE)
write.csv(lgg_methy, "survival-results/lgg_methy.csv", row.names = TRUE)
