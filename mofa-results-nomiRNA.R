setwd("~/Desktop/artigo")

source("setup.R")
### Load data
clinical <- read.csv("DataSets/2_analysis/clinical.csv", row.names = 1)

survival <- read.csv("DataSets/2_analysis/survival_data.csv", row.names = 1)

mutations <- read.csv("DataSets/2_analysis/mutations.csv", row.names = 1)
dna <- read.csv("DataSets/2_analysis/dna.csv", row.names = 1)
mrna <- read.csv("DataSets/2_analysis/mrna.csv", row.names = 1)



info.methy <- read.csv("DataSets/processed_assays/info_methylation.csv", row.names = 1)
info.mirna <- read.csv("DataSets/processed_assays/info_rna_mirna.csv", row.names = 1)


omics.list <- list(Mutations = mutations, Methylation = dna, mRNA = mrna)




## Run Model
source("run_mofa.R")

output <- run_mofa(omics.list, n_runs = 10, num_factors = NULL, 
                   likelihoods = c("bernoulli", "gaussian", "gaussian"), 
                   clinical_data = clinical) 



#### Save Results
model <- output$model_output 
print(model)
print(output$time_taken)





#### Read Z and W 
factors <- as.matrix(read.csv("mofa-output-nomiRNA/factors.csv", row.names = 1))

loadings <- list(
  Mutations = as.matrix(read.csv("mofa-output-nomiRNA/loadings_mutations.csv", row.names = 1)),
  Methylation = as.matrix(read.csv("mofa-output-nomiRNA/loadings_Methylation.csv", row.names = 1)),
  mRNA = as.matrix(read.csv("mofa-output-nomiRNA/loadings_mrna.csv", row.names = 1))
)

### Results shown in the paper 

## Variance Decomposition

p__variance_decomposition <- local({
  p1 <- plot_variance_explained(model, plot_total = T)[[1]]
  p1 <- p1$data
  
  p1$value <- round(p1$value, 1)
  colnames(p1) <- c("Factor", "View", "Value", "Group")
  p1$Factor <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", p1$Factor)
  p1$View <- factor(p1$View, levels = c("mRNA", "Methylation", "Mutations"))
  
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
ggsave(filename = "mofa-results-nomiRNA/variance-decomposition.pdf",
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



