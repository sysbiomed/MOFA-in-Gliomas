### Some setups
library(ggpubr)
library(survival)
library(ComplexHeatmap)
library(circlize)
library(ChAMP)
library(MOFA2)
library(survival)
library(survminer)
library(dplyr)
library(fgsea)
library(missMethyl)
library(clusterProfiler)
library(ReactomePA)
library(ReactomeContentService4R)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(DESeq2)
library(edgeR)
library(ggrepel)
library(ENmix)


text_size <- 2
my_theme <- theme(
  axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
  axis.text.x = element_text(size=rel(text_size), vjust=0.5, color='black'),
  axis.title.y=element_text(size=rel(text_size), hjust=1, color='black'),
  axis.title.x=element_text(size=rel(text_size), vjust=1, color='black'),
  legend.text = element_text(size = rel(text_size), color = 'black'),
  legend.title = element_text(size = rel(text_size), color = 'black'),
  panel.grid.minor = element_blank(), 
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  panel.background = element_blank(),
) 

