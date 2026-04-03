source("setup.R")


## Load data 
methy <- read.csv("DataSets/original_assays/assay_methylation_for_dge.csv", 
               row.names = 1)
methy_beta <- t(M2B(methy))


mrna <- t(as.matrix(read.csv("DataSets/original_assays/assay_rna_coding.csv", row.names = 1)))


classification <- as.data.frame(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", 
                                         row.names = 1))

labels <- local({
  classification$classification.2021 <- ifelse(classification $classification.2021 == "glioblastoma", "GBM", 
                                               ifelse(classification $classification.2021  == "astrocytoma", "ASTRO",
                                                      ifelse(classification $classification.2021  == "oligodendroglioma", "OLIGO",
                                                             "UNCLASS")))
  classification <- classification[,-c(1, 2), drop = FALSE]
  missing_samples <- setdiff(colnames(methy_beta), rownames(classification))
  new_rows <- data.frame(classification.2021 = rep("UNCLASS", length(missing_samples)), 
                         row.names = missing_samples)
  classification <- rbind(classification, new_rows)
  labels <- classification[colnames(methy_beta), "classification.2021"]
  labels
})


### DGE on dna
gbm_vs_astro_methy <- champ.DMP(beta = methy_beta,
                   pheno = labels,
                   compare.group = c("GBM", "ASTRO"),
                   adjPVal = 0.05 ,
                   adjust.method = "BH",
                   arraytype = "450K")
# -1 for ASTRO and 1 for GBM


gbm_vs_oligo_methy <- champ.DMP(beta = methy_beta,
                              pheno = labels,
                              compare.group = c("GBM", "OLIGO"),
                              adjPVal = 0.05 ,
                              adjust.method = "BH",
                              arraytype = "450K")
# -1 for GBM and 1 for OLIGO

astro_vs_oligo_methy <- champ.DMP(beta = methy_beta,
                                  pheno = labels,
                                  compare.group = c("ASTRO", "OLIGO"),
                                  adjPVal = 0.05 ,
                                  adjust.method = "BH",
                                  arraytype = "450K")
# -1 for ASTRO and 1 for OLIGO 


dge_2_save_methy <- function(df){

  df$diffexpressed <- "NO"
  df$diffexpressed[!is.na(df$logFC)  & df$logFC > 0.3 & df$P.Value < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -0.3 & df$P.Value < 0.05] <- "DOWN"

  return (df)
}

gbm_astro_df <- dge_2_save_methy(gbm_vs_astro_methy$ASTRO_to_GBM)

oligo_gbm_df <- dge_2_save_methy(gbm_vs_oligo_methy$GBM_to_OLIGO)
oligo_gbm_df$diffexpressed <- ifelse(
  oligo_gbm_df$diffexpressed == "UP",
  "DOWN",
  ifelse(
    oligo_gbm_df$diffexpressed == "DOWN",
    "UP",
    oligo_gbm_df$diffexpressed
  )
)

oligo_astro_df <- dge_2_save_methy(astro_vs_oligo_methy$ASTRO_to_OLIGO)
oligo_astro_df$diffexpressed <- ifelse(
  oligo_astro_df$diffexpressed == "UP",
  "DOWN",
  ifelse(
    oligo_astro_df$diffexpressed == "DOWN",
    "UP",
    oligo_astro_df$diffexpressed
  )
)

write.csv(gbm_astro_df, 
          "dge-results/gbm_vs_astro_methy.csv", row.names = TRUE)
write.csv(oligo_gbm_df, 
          "dge-results/oligo_vs_gbm_methy.csv", row.names = TRUE)
write.csv(oligo_astro_df,
          "dge-results/astro_vs_oligo_methy.csv", row.names = TRUE)



### DGE on mrna
# normalization as in https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# the difference for what we did is that we do not perform voom stabilization
keep <- rowSums(cpm(mrna) > 1) >= 5
mrna <- mrna[keep,]

mrna <- DGEList(counts = mrna)
mrna <- normLibSizes(mrna, method = "TMM")

mrna$samples$group <- classification[rownames(mrna$samples), "classification.2021"]


# design matrix
design <- model.matrix(~ 0 + group, data = mrna$samples)
colnames(design) <- levels(factor(mrna$samples$group))
design


# dispersion
mrna <- estimateDisp(mrna, design, robust=TRUE)
mrna$common.dispersion


# fit 
fit <- glmQLFit(mrna, design, robust=TRUE)
head(fit$coefficients)

# the contrasts we want to analyse
con1 <- makeContrasts(GBM - ASTRO, levels=design)
con2 <- makeContrasts(GBM- OLIGO, levels=design)
con3 <- makeContrasts(ASTRO - OLIGO, levels=design)

# do the test
qlf1 <- glmQLFTest(fit, contrast=con1)
qlf2 <- glmQLFTest(fit, contrast=con2)
qlf3 <- glmQLFTest(fit, contrast=con3)


dge_2_save_mrna <- function(df){
  
  df$diffexpressed <- "NO"
  df$diffexpressed[!is.na(df$logFC)  & df$logFC > 1 & df$PValue < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$PValue < 0.05] <- "DOWN"
  
  return (df)
}


## Save 
write.csv(dge_2_save_mrna(qlf1$table), "dge-results/gbm_vs_astro_mrna.csv", row.names = TRUE)
write.csv(dge_2_save_mrnaqlf2$table), "dge-results/gbm_vs_oligo_mrna.csv", row.names = TRUE)
write.csv(dge_2_save_mrnaqlf3$table), "dge-results/astro_vs_oligo_mrna.csv", row.names = TRUE)


